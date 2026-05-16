"""Tests for aggregation helper functions and timepoint auto-detection."""

import pytest

from edge_spot_analyser.aggregation import (
    extract_timestamp,
    extract_wellnumber,
    extract_xy,
)

# ---------------------------------------------------------------------------
# extract_xy
# ---------------------------------------------------------------------------


class TestExtractXY:
    def test_xy_format(self):
        assert extract_xy("WellD10_Channel405_Seq0009-MaxIP_XY9_405.tif") == 9

    def test_4digit_format(self):
        # 4-digit: raw position returned (0-indexed)
        assert (
            extract_xy("Plate000_WellB02_Channel405,561,488,640_Seq0000-MaxIP_B02_0007_405.tif")
            == 7
        )

    def test_4digit_zero(self):
        assert extract_xy("WellB02_Channel405,561,488,640_Seq0000-MaxIP_B02_0000_405.tif") == 0

    def test_named_channel_xy(self):
        assert extract_xy("WellB02_Channel20x GFP,20x DAPI_Seq0000-MaxIP_XY3_20x DAPI.tif") == 3

    def test_no_match_raises(self):
        with pytest.raises(ValueError):
            extract_xy("garbage_filename.tif")


# ---------------------------------------------------------------------------
# extract_timestamp
# ---------------------------------------------------------------------------


class TestExtractTimestamp:
    def test_lowercase_t1(self):
        assert (
            extract_timestamp("WellB02_Channel405,561,488,640_Seq0000-MaxIP_t1_T1_B02_0000_405.tif")
            == 1
        )

    def test_lowercase_t3(self):
        assert (
            extract_timestamp("WellB03_Channel405,561,488,640_Seq0006-MaxIP_t3_T1_B03_0002_488.tif")
            == 3
        )

    def test_no_timepoint_returns_zero(self):
        assert (
            extract_timestamp(
                "Plate000_WellB02_Channel405,561,488,640_Seq0000-MaxIP_B02_0000_405.tif"
            )
            == 0
        )

    def test_crop_t1_no_experimental_timepoint(self):
        """crop_T1 has uppercase T only — no experimental timepoint."""
        assert (
            extract_timestamp("WellE03_Channel405,561,488,640_Seq0003-MaxIP_crop_T1_XY8_405.tif")
            == 0
        )

    def test_named_channel_no_timepoint(self):
        assert (
            extract_timestamp("WellB02_Channel20x GFP,20x DAPI_Seq0000-MaxIP_XY1_20x DAPI.tif") == 0
        )


# ---------------------------------------------------------------------------
# extract_wellnumber
# ---------------------------------------------------------------------------


class TestExtractWellnumber:
    def test_standard(self):
        assert (
            extract_wellnumber("WellD10_Channel405,561,488,640_Seq0009-MaxIP_XY9_405.tif") == "D10"
        )

    def test_with_plate_prefix(self):
        assert (
            extract_wellnumber(
                "Plate000_WellB02_Channel405,561,488,640_Seq0000-MaxIP_B02_0000_405.tif"
            )
            == "B02"
        )

    def test_named_channel(self):
        assert (
            extract_wellnumber("WellG03_Channel20x GFP,20x DAPI_Seq0000-MaxIP_XY2_20x GFP.tif")
            == "G03"
        )

    def test_no_match_raises(self):
        with pytest.raises(ValueError):
            extract_wellnumber("garbage.tif")
