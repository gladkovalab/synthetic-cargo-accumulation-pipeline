"""Tests for figure aggregation timepoint-aware loading."""

import pandas as pd
import pytest

from edge_spot_analyser.figure_aggregation import FigureAggregator


class TestLoadWellData:
    """Test _load_well_data with and without timepoint."""

    @pytest.fixture()
    def results_dir(self, tmp_path):
        """Create a mock results directory with static CSV files."""
        date_dir = tmp_path / "221216"
        date_dir.mkdir()

        # Normal static file (all timepoints mixed — should not exist for
        # time-varying dates, but test fallback anyway)
        df_all = pd.DataFrame({"B02": [0.1, 0.2], "B03": [0.3, 0.4]}, index=[1, 2])
        df_all.index.name = "XY"
        df_all.to_csv(date_dir / "edge_spot_fraction_of_total_miro_static.csv")

        # Per-timepoint static files
        df_t1 = pd.DataFrame({"B02": [0.11, 0.21], "B03": [0.31, 0.41]}, index=[1, 2])
        df_t1.index.name = "XY"
        df_t1.to_csv(date_dir / "edge_spot_fraction_of_total_miro_static_t1.csv")

        df_t2 = pd.DataFrame({"B02": [0.12, 0.22], "B03": [0.32, 0.42]}, index=[1, 2])
        df_t2.index.name = "XY"
        df_t2.to_csv(date_dir / "edge_spot_fraction_of_total_miro_static_t2.csv")

        return tmp_path

    @pytest.fixture()
    def aggregator(self, results_dir, tmp_path):
        """Create a FigureAggregator pointing at the mock results."""
        # Need a minimal Excel config — create one with an empty figure sheet
        config_path = tmp_path / "config.xlsx"
        with pd.ExcelWriter(config_path) as writer:
            pd.DataFrame({"Date": [], "ctrl": []}).to_excel(writer, sheet_name="Fig1", index=False)
            pd.DataFrame({"Date": [], "Threshold correction factor": []}).to_excel(
                writer, sheet_name="Otsu_params", index=False
            )
        return FigureAggregator(
            results_dir=results_dir,
            config_path=config_path,
            output_dir=tmp_path / "figures",
        )

    def test_load_without_timepoint(self, aggregator):
        series = aggregator._load_well_data("221216", "B02", "edge_spot_fraction_of_total_miro_static.csv")
        assert series is not None
        assert list(series.values) == [0.1, 0.2]

    def test_load_with_timepoint_uses_tp_file(self, aggregator):
        series = aggregator._load_well_data(
            "221216", "B02", "edge_spot_fraction_of_total_miro_static.csv", timepoint="t1"
        )
        assert series is not None
        assert list(series.values) == [0.11, 0.21]

    def test_load_with_timepoint_t2(self, aggregator):
        series = aggregator._load_well_data(
            "221216", "B03", "edge_spot_fraction_of_total_miro_static.csv", timepoint="t2"
        )
        assert series is not None
        assert list(series.values) == [0.32, 0.42]

    def test_load_with_missing_timepoint_falls_back(self, aggregator):
        """If the per-timepoint file doesn't exist, fall back to the base file."""
        series = aggregator._load_well_data(
            "221216", "B02", "edge_spot_fraction_of_total_miro_static.csv", timepoint="t99"
        )
        assert series is not None
        assert list(series.values) == [0.1, 0.2]

    def test_load_missing_well_returns_none(self, aggregator):
        series = aggregator._load_well_data("221216", "Z99", "edge_spot_fraction_of_total_miro_static.csv")
        assert series is None

    def test_load_missing_date_returns_none(self, aggregator):
        series = aggregator._load_well_data("999999", "B02", "edge_spot_fraction_of_total_miro_static.csv")
        assert series is None
