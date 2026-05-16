# Synthetic cargo accumulation pipeline

Image-analysis pipeline that quantifies the spatial distribution of a Miro1-based synthetic cargo across fluorescence microscopy fields: nuclear segmentation, perinuclear-ring intensity (with Gini, moments, mass-displacement), and detection of "edge spots" outside the ring. Produces per-date CSV outputs and cross-date figure tables, plus an optional cell-level filter for an `_561` co-stain channel.

Paper companion to the figures in https://www.biorxiv.org/content/10.1101/2024.09.13.612963v1 — a Python re-implementation of the original CellProfiler pipeline used in the preprint (spec preserved in `pipeline_files/`). The Python port was undertaken during the revision process so we could re-run with different parameters and add additional metrics; see "Differences from the CellProfiler pipeline" below.

Sister repos under [`gladkovalab`](https://github.com/gladkovalab):

- [`micropattern-cell-analysis`](https://github.com/gladkovalab/micropattern-cell-analysis) — micropattern-constrained imaging analysis (wedge-r profiles, slab metrics)
- [`synthetic-cargo-particle-tracking`](https://github.com/gladkovalab/synthetic-cargo-particle-tracking) — single-particle tracking of the same cargo, TRAK isoform conditions


## Overview

For each date the pipeline:

- Loads matched `_405.tif` (Hoechst, nuclear), `_488.tif` (Miro1 synthetic cargo), and optional `_561.tif` (co-stain marker) files across multiple wells and XY positions.
- Segments nuclei in the Hoechst channel via 3-class Otsu thresholding (with a smoothing pre-step).
- Expands the nuclear mask by 10 px to define a perinuclear region for measurement, and by 15 px to define a separate exclusion mask.
- Detects "edge spots" — bright Miro1 puncta *outside* the perinuclear region — via Robust Background thresholding on the masked 488 channel.
- Measures perinuclear-region intensities (Gini coefficient, moments, mean intensity) and the edge-spot peripheral accumulation score (edge-spot intensity as a fraction of total Miro).
- Aggregates per-date and across-date results into figure-ready CSVs based on a condition-to-(date, well) mapping in an Excel workbook.

Optional 561 co-stain mode: when an `_561.tif` is present alongside `_405`/`_488`, per-cell perinuclear 561 intensities are added to `Perinuclear_region.csv`, and the `edge-spot-cell-filter` CLI can later drop non-control cells below a per-date control percentile.

### Differences from the CellProfiler pipeline used in the preprint

- A smoothing module (median filter) before the Otsu step.
- Border nuclei are kept for perinuclear masking (so spots near border nuclei are properly excluded) but removed for Gini analysis (where truncated perinuclear rings are problematic).
- An intensity-based peripheral accumulation score: total edge-spot intensity as a fraction of total MIRO intensity (the per-image score reported in the paper), replacing the original count/nucleus edge-spot fraction. (The legacy count/nucleus, per-nucleus and edge-to-perinuclear variants were computed in earlier revisions but dropped in v1.0.0 — only the paper metric is emitted.)
- Automatically generates the necessary figures from a spreadsheet mapping condition to (date, wellnumber).
- Optional 561-channel co-staining: when a `_561.tif` is present alongside the `_405`/`_488` pair, the pipeline measures per-cell perinuclear 561 intensities. `edge-spot-cell-filter` can then drop non-control cells whose 561 expression falls below a per-date control percentile threshold — useful when the synthetic-cargo experiment uses an expression marker and infection efficiency varies across replicates.


## Quick Start

### Prerequisites

Install uv if not already installed:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
# Or on macOS: brew install uv
```

### Installation

```bash
uv sync
```

### Run the Pipeline

```bash
# Process all dates (includes per-date aggregation)
uv run edge-spot-pipeline -i inputs/ -o results/

# Process specific date with parallel workers
uv run edge-spot-pipeline -i inputs/ -o results/ -d 241223 -w 10

# Use custom parameters from Excel config
uv run edge-spot-pipeline -i inputs/ -o results/ -c config.xlsx

# Restrict --figures-only to a subset of figure sheets (rest of the workbook is ignored)
uv run edge-spot-pipeline -i inputs/ -o results/ -c config.xlsx \
    --figures-only --figure-sheets FigS8A FigS8B

# Skip aggregation (only run image processing)
uv run edge-spot-pipeline -i inputs/ -o results/ --skip-aggregate
```

### Run Aggregation Separately

```bash
# Per-date aggregation (edge spot fractions, Gini, etc.)
uv run edge-spot-aggregate -i results/241223/

# Figure aggregation (combine across dates for publication figures).
# --figure-sheets is optional; default = every figure sheet in the workbook.
uv run edge-spot-figure-aggregate -r results/ -c config.xlsx \
    --figure-sheets FigS8A FigS8B
```

### Generate Plots

```bash
# Plot all normalized figure CSVs
uv run edge-spot-plot -i results/figures/

# Plot a single figure
uv run edge-spot-plot -i results/figures/Fig1D_edge_spot_normalized.csv
```

If `<Figure>_cell_counts.csv` is present in the figures directory (always
emitted by `edge-spot-figure-aggregate`), the corresponding `<Figure>_gini_*`
plot will be annotated with the per-condition cell count.

### Filter cells by 561-marker brightness (optional)

When the input data has a `_561.tif` co-stain, the pipeline records per-cell
perinuclear 561 metrics in `Perinuclear_region.csv`. To restrict the analysis
to cells whose 561 intensity exceeds a per-date control percentile, run:

```bash
uv run edge-spot-cell-filter \
    -r results/ -o results_filtered/ -c config.xlsx \
    --figure-sheets FigS8A FigS8B FigS8C FigS8D FigS8E \
    --control noT
# defaults: --filter-column Intensity_MeanIntensity_Marker561, --percentile 90

# Then re-aggregate and plot against the filtered mirror:
for d in results_filtered/*/; do uv run edge-spot-aggregate -i "$d"; done
uv run edge-spot-figure-aggregate -r results_filtered/ -c config.xlsx \
    --figure-sheets FigS8A FigS8B FigS8C FigS8D FigS8E
uv run edge-spot-plot -i results_filtered/figures/
```

Each date's threshold is the `--percentile` of the named `--control`
condition's per-cell `--filter-column` values for that date. Control cells are
kept unchanged; non-control cells from any of the listed `--figure-sheets`
wells are dropped if they fall below the threshold. Cells from wells that
don't appear in any of those sheets are passed through. The output is a
drop-in mirror of the input results dir (`Perinuclear_region.csv`,
`Nuclei.csv`, `Expand_Nuclei.csv`, `All_measurements.csv` are rewritten;
`edge_spots.csv` and `Image.csv` are copied through). A `filter_summary.csv`
at the top level records thresholds and per-well kept/total counts.

Two behaviours of the filter to be aware of:

- **Control wells are passed through unfiltered.** The control condition
  defines the noise-floor distribution used to compute the threshold, so all
  of its cells are kept. Only non-control cells in the listed figure sheets
  are thresholded; cells in wells absent from those sheets pass through.
- **Edge-spot detections are not re-thresholded per cell.** Edge spots are
  detected as standalone objects and are not tagged with a parent cell, so
  `edge_spots.csv` is copied through unchanged and the edge-spot columns in
  `All_measurements.csv` (`edge_spot_count`, `intensity_total`,
  `fraction_of_total_miro`) still reflect every spot in the image, while the
  `Nuclei` counts reflect the filtered cell count. This is an approximation
  that holds when low-561 cells contribute few or no spots.

## Input Data Format

The pipeline expects TIF files organized in date-based subdirectories:

```
inputs/
  ├── 241223/
  │   ├── WellB02_Channel405,561,488,640_Seq0000-MaxIP_XY1_405.tif
  │   ├── WellB02_Channel405,561,488,640_Seq0000-MaxIP_XY1_488.tif
  │   └── ...
  └── 241217/
      └── ...
```

**Filename format**: `Well{WELL}_Channel{...}_Seq{SEQ}[-MaxIP]_XY{XY}_{CHANNEL}.tif`
- `WELL`: Well identifier (e.g., B02, D10)
- `SEQ`: Sequence number (e.g., 0000, 0009)
- `-MaxIP`: Optional (maximum intensity projection indicator)
- `XY`: XY position number (1-9)
- `CHANNEL`: `405` for Hoechst (nuclear), `488` for MIRO160mer (mitochondrial), `561` (optional, used for the cell filter)

Files are automatically matched by `(well, sequence, XY)` tuple. When a
`_561.tif` is present for the same key it is loaded automatically and per-cell
perinuclear 561 metrics are added to `Perinuclear_region.csv`; missing 561
files are tolerated silently.

## Parameter Configuration

Parameters can be configured via JSON file or Excel spreadsheet.

### JSON Configuration

```json
{
  "default": {
    "otsu_correction_factor": 0.45,
    "diameter_min": 30,
    "diameter_max": 100,
    "min_distance": 30,
    "edge_spot_diameter_min": 3,
    "edge_spot_diameter_max": 80,
    "edge_spot_correction_factor": 2.0
  },
  "241223": {
    "otsu_correction_factor": 0.50
  }
}
```

### Excel Configuration

The pipeline can read parameters from an Excel file with an `Otsu_params` sheet containing per-date settings and exclusions.

### Key Parameters

**Nuclei Segmentation**
- **`otsu_correction_factor`** (default: 0.45): Multiplier for Otsu threshold. Lower = more nuclei detected.
- **`diameter_min`** / **`diameter_max`** (default: 30/100): Nucleus diameter range in pixels
- **`min_distance`** (default: 30): Minimum distance between nuclei for declumping

**Edge Spot Detection**
- **`edge_spot_correction_factor`** (default: 2.0): Multiplier for threshold. Higher = fewer spots detected.
- **`edge_spot_diameter_min`** / **`edge_spot_diameter_max`** (default: 3/80): Spot diameter range in pixels

## Output Files

### Per-Date Output (from pipeline)

```
results/
  └── 241223/
      ├── Nuclei.csv                    # Nuclear measurements (Hoechst intensity)
      ├── Expand_Nuclei.csv             # Expanded nuclei (MIRO160mer, Gini, moments)
      ├── Perinuclear_region.csv        # Perinuclear ring measurements
      ├── edge_spots.csv                # Edge spot measurements
      ├── Image.csv                     # Image-level statistics
      ├── All_measurements.csv          # Combined summary
      ├── edge_spot_fraction_of_total_miro_static.csv  # Peripheral accumulation score per well
      ├── GINI_Gini_MIRO160mer_fov_median_static.csv  # Gini per FOV
      └── ...
```

### Figure Output (from figure aggregation)

```
results/
  └── figures/
      ├── Fig1D_edge_spot_raw.csv        # Raw values for Figure 1D
      ├── Fig1D_edge_spot_normalized.csv # Normalized (per-date, control=1.0)
      ├── Fig1D_gini_raw.csv
      ├── Fig1D_gini_normalized.csv
      └── ...
```

### Key Measurements

- **Gini coefficient** (`GINI_Gini_MIRO160mer`): Inequality measure (0=uniform, 1=concentrated)
- **Mass displacement** (`Intensity_MassDisplacement_MIRO160mer`): Distance between intensity and geometric centroids
- **Moments**: Mean, standard deviation, skewness, kurtosis
- **Peripheral accumulation score** (`edge_spot_fraction_of_total_miro`): total edge-spot MIRO intensity / total MIRO intensity per field of view (normalized to control per replicate for the paper figures)

## Algorithm Details

### Nuclei Segmentation (Module 7)
- **Method**: 3-class Otsu thresholding (uses upper threshold)
- **Pre-smoothing**: median filter (3 px artifact diameter) applied to the Hoechst channel before thresholding (CellProfiler Smooth module)
- **Threshold smoothing**: additional Gaussian (σ=1.0, from threshold smoothing scale=2.0) applied internally before the Otsu step
- **Declumping**: Intensity-based watershed with min_distance=30px
- **Size filtering**: 30-100 pixel diameter
- **Border removal**: Yes

### Perinuclear Region (Modules 8, 9, 14)
- **Expand_Nuclei**: Nuclei expanded by 10 pixels (for measurements)
- **Expand_Nuclei_for_mask**: Nuclei expanded by 15 pixels (for masking only)
- **Perinuclear_region**: Ring between eroded nuclei (1px erosion) and 10px expansion

### Edge Spot Detection (Module 11 - Robust Background)
1. Remove outliers: lowest 30%, highest 10%
2. Calculate robust_mean = MEDIAN of remaining pixels
3. Calculate robust_std = STD of remaining pixels
4. Threshold = (robust_mean + 2×robust_std) × correction_factor
5. Clip to [0.0035, 1.0]
6. Apply Gaussian smoothing (σ=0.65) to image before thresholding
7. Simple connected component labeling (no declumping)

### Gini Coefficient (Module 19)
```python
flattened = sort(pixels)
npix = len(pixels)
normalization = abs(mean(pixels)) * npix * (npix - 1)
kernel = (2 * arange(1, npix+1) - npix - 1) * abs(pixels)
gini = sum(kernel) / normalization
```

## File Structure

```
.
├── src/edge_spot_analyser/
│   ├── pipeline.py            # Main pipeline (CLI: edge-spot-pipeline)
│   ├── segmentation.py        # Segmentation algorithms
│   ├── measurements.py        # Measurement functions
│   ├── io_utils.py            # File I/O, CSV export
│   ├── aggregation.py         # Per-date aggregation (CLI: edge-spot-aggregate)
│   ├── figure_aggregation.py  # Cross-date figure tables (CLI: edge-spot-figure-aggregate)
│   ├── cell_filter.py         # Per-cell 561 filter (CLI: edge-spot-cell-filter)
│   └── plotting.py            # Figure plotting (CLI: edge-spot-plot)
├── tests/                     # Unit tests
├── pipeline_files/            # Original CellProfiler pipeline JSON (reference)
├── config_V3.xlsx             # Round 1 config (per-date Otsu params, Aggregate, Fig* sheets)
├── config_round2_V3.xlsx      # Round 2 config
├── pyproject.toml             # Package + uv build config
├── uv.lock                    # Pinned dependency versions
└── README.md
```

`inputs/`, `results/`, and `plots/` are *user-created* working directories — not part of the repo. The pipeline reads the former and writes the latter two relative to wherever you invoke `uv run edge-spot-*` from.

## License

MIT License. See [LICENSE](LICENSE) for details.