"""Cell-level filter that thresholds on a per-control percentile.

For each date in a results directory, compute a threshold from the cells of a
designated control condition (e.g. ``noT``, ``noV``) on a chosen per-cell metric
(default: the perinuclear mean of the 561 marker channel) and drop every
non-control cell that falls below it. Cells in control wells are kept
unchanged. Cells whose well does not appear in any of the given figure sheets
are passed through.

Outputs a filtered mirror of the input results dir: filtered
``Perinuclear_region.csv`` / ``Nuclei.csv`` / ``Expand_Nuclei.csv`` per date,
plus an updated ``All_measurements.csv`` whose nuclei counts and
``edge_spot_intensity_per_nucleus`` reflect the filter. ``edge_spots.csv`` and
``Image.csv`` are copied through unchanged. A ``filter_summary.csv`` summary
of thresholds and per-well kept/total counts is written at the top level.

Run via:

    edge-spot-cell-filter -r results/ -o results_filtered/ -c config.xlsx \
        --figure-sheets FigS8A FigS8B --control noT
"""
from __future__ import annotations

import argparse
import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from edge_spot_analyser.io_utils import detect_figure_sheets

logger = logging.getLogger(__name__)

DEFAULT_FILTER_COLUMN = "Intensity_MeanIntensity_Marker561"
DEFAULT_PERCENTILE = 90.0
FILENAME_WELL_RE = re.compile(r"^(?:Plate\d+_)?Well([A-Z]\d+)_")


def _parse_well(filename: str) -> str | None:
    m = FILENAME_WELL_RE.match(str(filename))
    return m.group(1) if m else None


@dataclass
class CellFilterConfig:
    results_dir: Path
    output_dir: Path
    config_path: Path
    figure_sheets: list[str]
    control_condition: str
    filter_column: str = DEFAULT_FILTER_COLUMN
    percentile: float = DEFAULT_PERCENTILE


def _load_figure_mappings(
    config_path: Path, figure_sheets: list[str], control_condition: str
) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    """Build (long-form mapping, dict[date] -> set of control wells)."""
    excel = pd.ExcelFile(config_path)
    sheets = detect_figure_sheets(excel, restrict_to=figure_sheets)

    long_rows: list[dict] = []
    control_wells: dict[str, set[str]] = {}
    for sheet in sheets:
        df = pd.read_excel(config_path, sheet_name=sheet)
        date_col = df.columns[0]
        for _, r in df.iterrows():
            if pd.isna(r[date_col]):
                continue
            d = (
                str(int(r[date_col]))
                if not isinstance(r[date_col], str)
                else r[date_col].strip()
            )
            for cond in df.columns[1:]:
                v = r[cond]
                if pd.isna(v):
                    continue
                for w in str(v).split(","):
                    w = w.strip()
                    if not w:
                        continue
                    long_rows.append(
                        {"figure": sheet, "date": d, "condition": str(cond), "well": w}
                    )
                    if str(cond) == control_condition:
                        control_wells.setdefault(d, set()).add(w)
    return pd.DataFrame(long_rows), control_wells


def _filter_one_date(
    cfg: CellFilterConfig,
    date: str,
    figure_wells: set[str],
    control_wells: set[str],
) -> dict:
    src_dir = cfg.results_dir / date
    dst_dir = cfg.output_dir / date
    dst_dir.mkdir(parents=True, exist_ok=True)

    peri = pd.read_csv(src_dir / "Perinuclear_region.csv")
    if cfg.filter_column not in peri.columns:
        raise ValueError(
            f"{cfg.filter_column} not in Perinuclear_region.csv at {src_dir}. "
            f"Did the pipeline run with the 561 channel available?"
        )
    peri["_well"] = peri["FileName_Hoechst"].map(_parse_well)

    control_mask = peri["_well"].isin(control_wells)
    control_values = peri.loc[control_mask, cfg.filter_column].dropna().to_numpy()
    if control_values.size == 0:
        raise RuntimeError(
            f"No {cfg.control_condition} cells for date {date} "
            f"(expected wells: {sorted(control_wells)})."
        )
    threshold = float(np.percentile(control_values, cfg.percentile))

    in_figs = peri["_well"].isin(figure_wells)
    is_ctrl = peri["_well"].isin(control_wells)
    keep_mask = is_ctrl | (~in_figs) | (
        in_figs & ~is_ctrl & (peri[cfg.filter_column] >= threshold)
    )
    peri_kept = peri.loc[keep_mask].copy()
    keep_keys = set(zip(peri_kept["ImageNumber"], peri_kept["ObjectNumber"], strict=True))

    peri_kept.drop(columns=["_well"]).to_csv(
        dst_dir / "Perinuclear_region.csv", index=False
    )

    for name in ("Expand_Nuclei.csv", "Nuclei.csv"):
        src_csv = src_dir / name
        if not src_csv.exists():
            continue
        df = pd.read_csv(src_csv)
        keep_idx = df.apply(
            lambda r: (r["ImageNumber"], r["ObjectNumber"]) in keep_keys, axis=1
        )
        df.loc[keep_idx].to_csv(dst_dir / name, index=False)

    for name in ("edge_spots.csv", "Image.csv"):
        src_csv = src_dir / name
        if src_csv.exists():
            pd.read_csv(src_csv).to_csv(dst_dir / name, index=False)

    all_csv = src_dir / "All_measurements.csv"
    if all_csv.exists():
        all_df = pd.read_csv(all_csv, header=[0, 1])
        counts_by_image = peri_kept.groupby("ImageNumber").size().to_dict()
        new_n: list[int] = []
        new_per_nucleus: list[float] = []
        for _, row in all_df.iterrows():
            img_no = int(row[("Image", "ImageNumber")])
            n = int(counts_by_image.get(img_no, 0))
            new_n.append(n)
            intensity_total = float(row[("edge_spots", "intensity_total")])
            new_per_nucleus.append(intensity_total / n if n > 0 else 0.0)
        all_df[("Nuclei", "Number_Object_Number")] = new_n
        all_df[("Nuclei", "Count_Interior")] = new_n
        all_df[("edge_spots", "intensity_per_nucleus")] = new_per_nucleus
        all_df.to_csv(dst_dir / "All_measurements.csv", index=False)

    summary: dict = {"date": date, "threshold": threshold}
    for well in sorted(figure_wells):
        sub = peri[peri["_well"] == well]
        kept = peri_kept[peri_kept["_well"] == well]
        summary[f"n_total_{well}"] = int(len(sub))
        summary[f"n_kept_{well}"] = int(len(kept))
    return summary


def apply_cell_filter(cfg: CellFilterConfig) -> pd.DataFrame:
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    mapping, control_wells = _load_figure_mappings(
        cfg.config_path, cfg.figure_sheets, cfg.control_condition
    )
    if mapping.empty:
        raise RuntimeError(
            f"No (date, well) entries from figure sheets {cfg.figure_sheets} in {cfg.config_path}"
        )

    summaries: list[dict] = []
    for date in sorted(mapping["date"].unique()):
        figure_wells = set(mapping.loc[mapping["date"] == date, "well"])
        ctrl_for_date = control_wells.get(date, set())
        if not ctrl_for_date:
            logger.warning(
                "Date %s has no %s wells across figure sheets %s — skipping",
                date,
                cfg.control_condition,
                cfg.figure_sheets,
            )
            continue
        date_dir = cfg.results_dir / date
        if not date_dir.exists():
            logger.warning("Skipping date %s (no results dir at %s)", date, date_dir)
            continue
        s = _filter_one_date(cfg, date, figure_wells, ctrl_for_date)
        summaries.append(s)
        keep_pcts = []
        for well in sorted(figure_wells):
            tot = s.get(f"n_total_{well}", 0)
            kpt = s.get(f"n_kept_{well}", 0)
            tag = cfg.control_condition if well in ctrl_for_date else ""
            keep_pcts.append(
                f"{well}{tag}: {kpt}/{tot}"
                + (f" ({100.0 * kpt / tot:.1f}%)" if tot else "")
            )
        logger.info("%s: thr=%.5f | %s", date, s["threshold"], " | ".join(keep_pcts))

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(cfg.output_dir / "filter_summary.csv", index=False)
    return summary_df


def main() -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    p = argparse.ArgumentParser(
        description=(
            "Filter pipeline outputs by per-date control-condition percentile threshold. "
            "Default: drop cells whose perinuclear 561 mean intensity falls below the "
            "90th percentile of the control wells' perinuclear 561 mean."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  edge-spot-cell-filter -r results/ -o results_filtered/ -c config.xlsx \\
      --figure-sheets FigS8A FigS8B FigS8C FigS8D FigS8E --control noT

  edge-spot-cell-filter -r results_s14/ -o results_s14_filtered/ -c config.xlsx \\
      --figure-sheets FigS14A FigS14B FigS14C FigS14F FigS14G --control noV
""",
    )
    p.add_argument("--results", "-r", type=Path, required=True)
    p.add_argument("--output", "-o", type=Path, required=True)
    p.add_argument("--config", "-c", type=Path, required=True)
    p.add_argument(
        "--figure-sheets",
        nargs="+",
        required=True,
        help="Figure sheet names whose (date, well, condition) mapping defines the filter.",
    )
    p.add_argument(
        "--control",
        required=True,
        help="Condition column name treated as the negative control (e.g. noT, noV).",
    )
    p.add_argument(
        "--filter-column",
        default=DEFAULT_FILTER_COLUMN,
        help=f"Per-cell column to threshold on (default: {DEFAULT_FILTER_COLUMN})",
    )
    p.add_argument(
        "--percentile",
        type=float,
        default=DEFAULT_PERCENTILE,
        help=f"Percentile of control distribution used as threshold (default: {DEFAULT_PERCENTILE})",
    )

    args = p.parse_args()
    cfg = CellFilterConfig(
        results_dir=args.results,
        output_dir=args.output,
        config_path=args.config,
        figure_sheets=args.figure_sheets,
        control_condition=args.control,
        filter_column=args.filter_column,
        percentile=args.percentile,
    )
    apply_cell_filter(cfg)


if __name__ == "__main__":
    main()
