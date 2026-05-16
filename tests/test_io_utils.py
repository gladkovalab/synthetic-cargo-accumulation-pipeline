"""Tests for I/O utilities: filename parsing, pairing, config loading, and filtering."""

from pathlib import Path

from edge_spot_analyser.io_utils import (
    FileDiscovery,
    ImagePair,
    _parse_well_timepoint,
    parse_exclusions,
    should_exclude_image,
    should_include_image,
)

# ---------------------------------------------------------------------------
# _parse_filename
# ---------------------------------------------------------------------------


class TestParseFilename:
    """Test all filename patterns handled by FileDiscovery._parse_filename."""

    # Pattern 1: standard XY format
    def test_standard_xy(self):
        result = FileDiscovery._parse_filename(
            "WellD10_Channel405,561,488,640_Seq0009-MaxIP_XY9_405"
        )
        assert result == {
            "well": "D10",
            "sequence": "Seq0009",
            "xy": "XY9",
            "channel": "405",
        }

    def test_standard_xy_with_plate_prefix(self):
        result = FileDiscovery._parse_filename(
            "Plate000_WellB02_Channel405,561,488,640_Seq0000-MaxIP_XY3_488"
        )
        assert result == {
            "well": "B02",
            "sequence": "Seq0000",
            "xy": "XY3",
            "channel": "488",
        }

    def test_standard_xy_no_maxip(self):
        result = FileDiscovery._parse_filename("WellC04_Channel405,488_Seq0002_XY5_405")
        assert result == {
            "well": "C04",
            "sequence": "Seq0002",
            "xy": "XY5",
            "channel": "405",
        }

    def test_named_channel_xy(self):
        """Round 3 format: named channels like '20x DAPI' instead of '405'."""
        result = FileDiscovery._parse_filename(
            "WellB02_Channel20x GFP,20x DAPI_Seq0000-MaxIP_XY1_20x DAPI"
        )
        assert result == {
            "well": "B02",
            "sequence": "Seq0000",
            "xy": "XY1",
            "channel": "20x DAPI",
        }

    def test_named_channel_gfp(self):
        result = FileDiscovery._parse_filename(
            "WellB02_Channel20x GFP,20x DAPI_Seq0000-MaxIP_XY1_20x GFP"
        )
        assert result["channel"] == "20x GFP"

    # Pattern 2: crop_T1 / t prefix + XY
    def test_crop_t1_xy(self):
        result = FileDiscovery._parse_filename(
            "WellE03_Channel405,561,488,640_Seq0003-MaxIP_crop_T1_XY8_405"
        )
        assert result == {
            "well": "E03",
            "sequence": "Seq0003",
            "timepoint": None,  # crop_T1 has no experimental timepoint
            "xy": "XY8",
            "channel": "405",
        }

    def test_t1_xy_with_timepoint(self):
        result = FileDiscovery._parse_filename(
            "WellB03_Channel405,561,488,640_Seq0006-MaxIP_t2_T1_XY4_488"
        )
        assert result["timepoint"] == "t2"
        assert result["well"] == "B03"
        assert result["xy"] == "XY4"

    def test_t_prefix_without_uppercase_T(self):
        """Bare t prefix without _T1 suffix should still parse."""
        result = FileDiscovery._parse_filename(
            "WellB03_Channel405,561,488,640_Seq0006-MaxIP_t1_XY4_488"
        )
        assert result is not None
        assert result["timepoint"] == "t1"

    # Pattern 3: 4-digit position format
    def test_4digit_position(self):
        result = FileDiscovery._parse_filename(
            "WellF02_Channel405,561,488,640_Seq0001-MaxIP_F02_0007_405"
        )
        assert result == {
            "well": "F02",
            "sequence": "Seq0001",
            "xy": "XY8",  # 0007 -> (7 % 10) + 1 = 8
            "channel": "405",
        }

    def test_4digit_with_plate_prefix(self):
        result = FileDiscovery._parse_filename(
            "Plate000_WellB02_Channel405,561,488,640_Seq0000-MaxIP_B02_0000_405"
        )
        assert result["well"] == "B02"
        assert result["xy"] == "XY1"  # 0000 -> (0 % 10) + 1 = 1

    def test_4digit_named_channel(self):
        result = FileDiscovery._parse_filename(
            "WellB02_Channel20x GFP,20x DAPI_Seq0000-MaxIP_B02_0000_20x DAPI"
        )
        assert result["well"] == "B02"
        assert result["channel"] == "20x DAPI"
        assert result["xy"] == "XY1"

    # Pattern 4: 4-digit with crop_T1 / t prefix
    def test_4digit_crop_t1(self):
        result = FileDiscovery._parse_filename(
            "WellB02_Channel405,561,488,640_Seq0000-MaxIP_crop_T1_B02_0000_405"
        )
        assert result["timepoint"] is None
        assert result["well"] == "B02"

    def test_4digit_t1_timepoint(self):
        result = FileDiscovery._parse_filename(
            "WellB02_Channel405,561,488,640_Seq0000-MaxIP_t1_T1_B02_0000_405"
        )
        assert result["timepoint"] == "t1"
        assert result["well"] == "B02"
        assert result["xy"] == "XY1"

    def test_4digit_t3_timepoint(self):
        result = FileDiscovery._parse_filename(
            "WellB03_Channel405,561,488,640_Seq0006-MaxIP_t3_T1_B03_0002_488"
        )
        assert result["timepoint"] == "t3"
        assert result["well"] == "B03"
        assert result["xy"] == "XY3"  # 0002 -> (2 % 10) + 1 = 3

    def test_4digit_t_without_uppercase_T(self):
        result = FileDiscovery._parse_filename(
            "WellB02_Channel405,561,488,640_Seq0000-MaxIP_t2_B02_0001_405"
        )
        assert result is not None
        assert result["timepoint"] == "t2"
        assert result["xy"] == "XY2"

    # Edge cases
    def test_composite_image_no_channel_suffix_returns_none(self):
        """Composite TIFs (no channel suffix) should not parse."""
        result = FileDiscovery._parse_filename(
            "Plate000_WellB02_Channel405,561,488,640_Seq0000-MaxIP_B02_0000"
        )
        assert result is None

    def test_garbage_returns_none(self):
        assert FileDiscovery._parse_filename("not_a_valid_filename") is None

    def test_duplicate_channel_numbers(self):
        """Some plates have duplicate channels like 405,561,488,561."""
        result = FileDiscovery._parse_filename(
            "WellB04_Channel405,561,488,561_Seq0002-MaxIP_B04_0000_405"
        )
        assert result is not None
        assert result["well"] == "B04"


# ---------------------------------------------------------------------------
# _match_image_pairs  (timepoint collision)
# ---------------------------------------------------------------------------


class TestMatchImagePairs:
    """Verify that images with different timepoints don't collide."""

    def _make_paths(self, names: list[str], tmp_path: Path) -> list[Path]:
        paths = []
        for name in names:
            p = tmp_path / f"{name}.tif"
            p.touch()
            paths.append(p)
        return paths

    def test_different_timepoints_produce_separate_pairs(self, tmp_path):
        hoechst = self._make_paths(
            [
                "WellB02_Channel405,561,488,640_Seq0000-MaxIP_t1_T1_B02_0000_405",
                "WellB02_Channel405,561,488,640_Seq0000-MaxIP_t2_T1_B02_0000_405",
                "WellB02_Channel405,561,488,640_Seq0000-MaxIP_t3_T1_B02_0000_405",
            ],
            tmp_path,
        )
        miro = self._make_paths(
            [
                "WellB02_Channel405,561,488,640_Seq0000-MaxIP_t1_T1_B02_0000_488",
                "WellB02_Channel405,561,488,640_Seq0000-MaxIP_t2_T1_B02_0000_488",
                "WellB02_Channel405,561,488,640_Seq0000-MaxIP_t3_T1_B02_0000_488",
            ],
            tmp_path,
        )
        pairs = FileDiscovery._match_image_pairs(hoechst, miro, "221216")
        assert len(pairs) == 3
        timepoints = sorted(p.timepoint for p in pairs)
        assert timepoints == ["t1", "t2", "t3"]

    def test_no_timepoint_pairs_still_work(self, tmp_path):
        hoechst = self._make_paths(
            [
                "WellB02_Channel405,561,488,640_Seq0000-MaxIP_XY1_405",
                "WellB02_Channel405,561,488,640_Seq0000-MaxIP_XY2_405",
            ],
            tmp_path,
        )
        miro = self._make_paths(
            [
                "WellB02_Channel405,561,488,640_Seq0000-MaxIP_XY1_488",
                "WellB02_Channel405,561,488,640_Seq0000-MaxIP_XY2_488",
            ],
            tmp_path,
        )
        pairs = FileDiscovery._match_image_pairs(hoechst, miro, "220709")
        assert len(pairs) == 2
        assert all(p.timepoint is None for p in pairs)

    def test_unmatched_miro_ignored(self, tmp_path):
        hoechst = self._make_paths(["WellB02_Channel405,488_Seq0000-MaxIP_XY1_405"], tmp_path)
        miro = self._make_paths(
            [
                "WellB02_Channel405,488_Seq0000-MaxIP_XY1_488",
                "WellB02_Channel405,488_Seq0000-MaxIP_XY9_488",  # no matching hoechst
            ],
            tmp_path,
        )
        pairs = FileDiscovery._match_image_pairs(hoechst, miro, "test")
        assert len(pairs) == 1


# ---------------------------------------------------------------------------
# ImagePair
# ---------------------------------------------------------------------------


class TestImagePair:
    def test_image_id_without_timepoint(self):
        pair = ImagePair(
            hoechst_path=Path("h.tif"),
            miro_path=Path("m.tif"),
            well_number="B02",
            sequence="Seq0000",
            xy_position="XY1",
            date="220709",
        )
        assert pair.image_id == "220709_B02_Seq0000_XY1"

    def test_image_id_with_timepoint(self):
        pair = ImagePair(
            hoechst_path=Path("h.tif"),
            miro_path=Path("m.tif"),
            well_number="B02",
            sequence="Seq0000",
            xy_position="XY1",
            date="221216",
            timepoint="t2",
        )
        assert pair.image_id == "221216_B02_t2_Seq0000_XY1"

    def test_metadata_dict_includes_timepoint(self):
        pair = ImagePair(
            hoechst_path=Path("h.tif"),
            miro_path=Path("m.tif"),
            well_number="B03",
            sequence="Seq0006",
            xy_position="XY4",
            date="221216",
            timepoint="t1",
        )
        meta = pair.to_metadata_dict(1)
        assert meta["Metadata_Timepoint"] == "t1"
        assert meta["Metadata_Well"] == "B03"

    def test_metadata_dict_omits_timepoint_when_none(self):
        pair = ImagePair(
            hoechst_path=Path("h.tif"),
            miro_path=Path("m.tif"),
            well_number="B02",
            sequence="Seq0000",
            xy_position="XY1",
            date="220709",
        )
        meta = pair.to_metadata_dict(1)
        assert "Metadata_Timepoint" not in meta


# ---------------------------------------------------------------------------
# parse_exclusions
# ---------------------------------------------------------------------------


class TestParseExclusions:
    def test_xy_format(self):
        result = parse_exclusions("E06,XY9 & D05, XY6")
        assert result == {("E06", 9), ("D05", 6)}

    def test_4digit_format(self):
        result = parse_exclusions("C09,0007 & D04,0004")
        assert result == {("C09", 8), ("D04", 5)}  # 0-indexed -> 1-indexed

    def test_mixed_formats(self):
        result = parse_exclusions("E06,XY9 & C09,0007")
        assert ("E06", 9) in result
        assert ("C09", 8) in result

    def test_empty_string(self):
        assert parse_exclusions("") == set()

    def test_nan_string(self):
        assert parse_exclusions("nan") == set()

    def test_single_exclusion(self):
        result = parse_exclusions("F04, XY1")
        assert result == {("F04", 1)}


class TestShouldExcludeImage:
    def _make_pair(self, well: str, xy: str) -> ImagePair:
        return ImagePair(
            hoechst_path=Path("h.tif"),
            miro_path=Path("m.tif"),
            well_number=well,
            sequence="Seq0000",
            xy_position=xy,
            date="220709",
        )

    def test_excluded_pair(self):
        exclusions = {("E06", 9)}
        assert should_exclude_image(self._make_pair("E06", "XY9"), exclusions)

    def test_non_excluded_pair(self):
        exclusions = {("E06", 9)}
        assert not should_exclude_image(self._make_pair("B02", "XY1"), exclusions)

    def test_empty_exclusions(self):
        assert not should_exclude_image(self._make_pair("B02", "XY1"), set())


# ---------------------------------------------------------------------------
# _parse_well_timepoint
# ---------------------------------------------------------------------------


class TestParseWellTimepoint:
    def test_plain_well(self):
        assert _parse_well_timepoint("B03") == ("B03", None)

    def test_well_with_timepoint(self):
        assert _parse_well_timepoint("B03 t1") == ("B03", "t1")

    def test_well_with_higher_timepoint(self):
        assert _parse_well_timepoint("D04 t12") == ("D04", "t12")

    def test_whitespace_stripped(self):
        assert _parse_well_timepoint("  B03 t1  ") == ("B03", "t1")

    def test_non_timepoint_suffix_kept_as_well(self):
        """Suffixes that don't match t{N} are not treated as timepoints."""
        well, tp = _parse_well_timepoint("B03 extra")
        assert tp is None


# ---------------------------------------------------------------------------
# should_include_image  (timepoint-aware)
# ---------------------------------------------------------------------------


class TestShouldIncludeImage:
    def _make_pair(self, date: str, well: str, timepoint: str | None = None) -> ImagePair:
        return ImagePair(
            hoechst_path=Path("h.tif"),
            miro_path=Path("m.tif"),
            well_number=well,
            sequence="Seq0000",
            xy_position="XY1",
            date=date,
            timepoint=timepoint,
        )

    def test_none_inclusions_includes_everything(self):
        assert should_include_image(self._make_pair("220709", "B02"), None)

    def test_exact_match_no_timepoint(self):
        inclusions = {("220709", "B02", None)}
        assert should_include_image(self._make_pair("220709", "B02"), inclusions)

    def test_no_match(self):
        inclusions = {("220709", "B02", None)}
        assert not should_include_image(self._make_pair("220709", "C03"), inclusions)

    def test_exact_match_with_timepoint(self):
        inclusions = {("221216", "B03", "t1")}
        assert should_include_image(self._make_pair("221216", "B03", "t1"), inclusions)

    def test_wrong_timepoint_rejected(self):
        inclusions = {("221216", "B03", "t1")}
        assert not should_include_image(self._make_pair("221216", "B03", "t2"), inclusions)

    def test_none_timepoint_inclusion_matches_any_image_timepoint(self):
        """An inclusion with timepoint=None should match images regardless of their timepoint."""
        inclusions = {("221216", "B03", None)}
        assert should_include_image(self._make_pair("221216", "B03", "t2"), inclusions)
        assert should_include_image(self._make_pair("221216", "B03", None), inclusions)

    def test_timepoint_inclusion_does_not_match_no_timepoint_image(self):
        """An inclusion requiring t1 should not match an image without any timepoint."""
        inclusions = {("221216", "B03", "t1")}
        assert not should_include_image(self._make_pair("221216", "B03", None), inclusions)
