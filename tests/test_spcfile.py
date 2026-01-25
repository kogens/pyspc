from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from spcfile import SPCFile, SPCSubfile
from spcfile.spcfile import FLAG_EXPLICIT_X
from tests.fixtures_synthetic import choose_data_dir


@pytest.fixture(scope="session")
def data_dir(tmp_path_factory) -> Path:
    real = Path(__file__).parent / "data"
    return choose_data_dir(tmp_path_factory=tmp_path_factory, real_data_dir=real)


class TestSPCFileConstruction:
    """Test file loading and basic construction."""

    @pytest.mark.parametrize("filename", ["s_evenx.spc", "s_xy.spc", "ft-ir.spc"])
    def test_load_single_spectrum(self, data_dir: Path, filename: str) -> None:
        """Load single-spectrum files."""
        spc = SPCFile(data_dir / filename)
        assert isinstance(spc, SPCFile)
        assert spc.path == data_dir / filename

    @pytest.mark.parametrize("filename", ["m_evenz.spc", "nir.spc"])
    def test_load_multifile(self, data_dir: Path, filename: str) -> None:
        """Load multifile spectra."""
        spc = SPCFile(data_dir / filename)
        assert isinstance(spc, SPCFile)
        assert len(spc) > 1

    def test_4d_file_minimal_support(self, data_dir: Path) -> None:
        """4D W-plane files should load and expose minimal W metadata."""
        spc = SPCFile(data_dir / "4d_map.spc")
        assert spc.w_planes > 0
        # Per-subfile W coordinate
        assert spc.w.shape == (len(spc),)
        assert np.all(np.isfinite(spc.w))

        # Known basic properties for this sample file.
        assert spc.w_planes == 11
        assert len(spc) == 121

    def test_nonexistent_file(self) -> None:
        """Raise FileNotFoundError for missing files."""
        with pytest.raises(FileNotFoundError):
            SPCFile("does_not_exist.spc")

    def test_unsupported_version(self, data_dir: Path) -> None:
        """Reject unsupported file versions with clear error."""
        with pytest.raises(ValueError, match="Unsupported SPC version"):
            SPCFile(data_dir / "m_ordz.spc")  # 0x4D old format


class TestSharedXMode:
    """Test files with a single shared X axis (current implementation)."""

    @pytest.fixture(scope="class")
    def single_spc(self, data_dir: Path) -> SPCFile:
        """Load s_evenx.spc once for the entire test class."""
        return SPCFile(data_dir / "s_evenx.spc")

    @pytest.fixture(scope="class")
    def multi_spc(self, data_dir: Path) -> SPCFile:
        """Load m_evenz.spc once for the entire test class."""
        return SPCFile(data_dir / "m_evenz.spc")

    @pytest.mark.parametrize("filename", ["s_evenx.spc", "s_xy.spc", "m_evenz.spc", "ft-ir.spc"])
    def test_shared_x_mode_contract(self, data_dir: Path, filename: str) -> None:
        """Shared-X files expose a stable x/y array API."""
        spc = SPCFile(data_dir / filename)

        assert spc.has_shared_x is True

        assert spc.x.ndim == 1
        assert spc.y.ndim in (1, 2)
        assert spc.x.shape[0] == spc.y.shape[0]

        assert spc.x.dtype == np.float64
        assert spc.y.dtype == np.float64

        assert np.all(np.isfinite(spc.x))
        assert np.all(np.isfinite(spc.y))

        dx = np.diff(spc.x)
        assert np.all(dx > 0) or np.all(dx < 0)

        if spc.y.ndim == 2:
            assert spc.y.shape[1] == len(spc)

    def test_expected_shapes_for_representative_files(self, single_spc: SPCFile, multi_spc: SPCFile) -> None:
        """A couple of concrete files act as regression tests for shapes."""
        assert single_spc.x.shape == (1844,)
        assert single_spc.y.shape == (1844,)

        assert multi_spc.x.shape == (171,)
        assert multi_spc.y.shape == (171, 32)


class TestXModeBehavior:
    """Test explicit vs implicit X handling."""

    def test_implicit_evenly_spaced_x(self, data_dir: Path) -> None:
        """Implicit X files rely on the header endpoints/linspace."""
        spc = SPCFile(data_dir / "s_evenx.spc")
        assert not (spc.header["flags"] & FLAG_EXPLICIT_X)

        expected = np.linspace(spc.header["first_x"], spc.header["last_x"], spc.header["n_points"])
        assert np.allclose(spc.x, expected, rtol=1e-6, atol=0.0)

    def test_explicit_x_overrides_linspace(self, data_dir: Path) -> None:
        """Explicit global X arrays should differ from implied linspace."""
        spc = SPCFile(data_dir / "s_xy.spc")
        assert spc.header["flags"] & FLAG_EXPLICIT_X

        expected = np.linspace(spc.header["first_x"], spc.header["last_x"], spc.header["n_points"])
        assert not np.allclose(spc.x, expected, rtol=1e-6, atol=0.0)


class TestPerSubfileXYMode:
    """Test TXYXYS files where each subfile has its own X/Y arrays."""

    def test_m_xyxy_contract(self, data_dir: Path) -> None:
        spc = SPCFile(data_dir / "m_xyxy.spc")
        assert spc.has_shared_x is False
        assert spc.per_subfile_xy is True
        assert spc.explicit_x is True

        with pytest.raises(ValueError):
            _ = spc.x
        with pytest.raises(ValueError):
            _ = spc.y

        assert len(spc) == int(spc.header["n_subfiles"])

        sub0 = spc[0]
        assert isinstance(sub0, SPCSubfile)
        assert sub0.x.dtype == np.float64
        assert sub0.y.dtype == np.float64
        assert sub0.x.ndim == 1
        assert sub0.y.ndim == 1
        assert sub0.x.shape == sub0.y.shape
        assert sub0.x.shape[0] == int(sub0.subheader["n_points"])
        assert np.all(np.isfinite(sub0.x))
        assert np.all(np.isfinite(sub0.y))

    def test_single_subfile_xyxy_allows_x_y_access(self, data_dir: Path) -> None:
        """Single-subfile files with TXYXYS set should still allow spc.x/spc.y for convenience."""
        spc = SPCFile(data_dir / "ms.spc")
        assert len(spc) == 1
        assert spc.per_subfile_xy is True

        assert spc.x.ndim == 1
        assert spc.y.ndim == 1
        assert spc.x.shape == spc[0].x.shape
        assert spc.y.shape == spc[0].y.shape
        assert np.array_equal(spc.x, spc[0].x)
        assert np.array_equal(spc.y, spc[0].y)


class TestSPCFileIndexing:
    """Test indexing and iteration over subfiles."""

    @pytest.fixture(scope="class")
    def multi_spc(self, data_dir: Path) -> SPCFile:
        """Load m_evenz.spc once for the entire test class."""
        return SPCFile(data_dir / "m_evenz.spc")

    def test_len_single_spectrum(self, data_dir: Path) -> None:
        """Single spectrum files should have length 1."""
        spc = SPCFile(data_dir / "s_evenx.spc")
        assert len(spc) == 1

    def test_indexing_and_iteration_contract(self, multi_spc: SPCFile) -> None:
        """Indexing/iteration always yields SPCSubfile views consistent with x/y arrays."""
        assert len(multi_spc) == 32

        sub0 = multi_spc[0]
        assert isinstance(sub0, SPCSubfile)
        assert sub0.x.shape == multi_spc.x.shape
        assert sub0.y.shape == (multi_spc.x.shape[0],)

        # Check a few representative indices (start, middle, end).
        for i in (0, 5, len(multi_spc) - 1):
            assert np.array_equal(multi_spc[i].x, multi_spc.x)
            assert np.array_equal(multi_spc[i].y, multi_spc.y[:, i])
            assert isinstance(multi_spc[i].z, (int, float))

        subfiles = list(multi_spc)
        assert len(subfiles) == len(multi_spc)
        for i, sub in enumerate(subfiles):
            assert isinstance(sub, SPCSubfile)
            assert np.array_equal(sub.x, multi_spc[i].x)
            assert np.array_equal(sub.y, multi_spc[i].y)


class TestSPCSubfile:
    """Test SPCSubfile properties."""

    @pytest.fixture(scope="class")
    def multi_spc(self, data_dir: Path) -> SPCFile:
        """Load m_evenz.spc once for the entire test class."""
        return SPCFile(data_dir / "m_evenz.spc")

    def test_subfile_basic_contract(self, multi_spc: SPCFile) -> None:
        """A returned SPCSubfile has aligned x/y arrays and basic metadata."""
        sub = multi_spc[5]
        assert sub.x.ndim == 1
        assert sub.y.ndim == 1
        assert sub.x.shape[0] == sub.y.shape[0]
        assert np.all(np.isfinite(sub.x))
        assert np.all(np.isfinite(sub.y))
        assert isinstance(sub.z, (int, float))


class TestLogText:
    """Test log text parsing."""

    @pytest.fixture(scope="class")
    def ftir_spc(self, data_dir: Path) -> SPCFile:
        """Load ft-ir.spc once for the entire test class."""
        return SPCFile(data_dir / "ft-ir.spc")

    def test_log_text_present_and_readable(self, ftir_spc: SPCFile) -> None:
        """Log text is present and looks like human-readable key/value pairs."""
        assert ftir_spc.log is not None
        assert isinstance(ftir_spc.log, str)
        assert len(ftir_spc.log) > 0
        assert ("MODEL" in ftir_spc.log) or ("SCANS" in ftir_spc.log)
        assert ("\r\n" in ftir_spc.log) or ("\n" in ftir_spc.log)

    def test_log_header_skipped(self, ftir_spc: SPCFile) -> None:
        """Log text should start with human-readable records, not binary header."""
        assert ftir_spc.log is not None
        text = ftir_spc.log
        assert text[0] != "\x00"
        assert text.startswith("MODEL")
        assert len(text) == 376  # Exact log length for this file

    def test_log_absent_or_none(self, data_dir: Path) -> None:
        """Files without log should have None."""
        spc = SPCFile(data_dir / "s_evenx.spc")
        assert spc.log is None or isinstance(spc.log, str)


class TestKnownFileProperties:
    """Test specific files with known golden values."""

    def test_s_evenx_properties(self, data_dir: Path) -> None:
        """Test known properties of s_evenx.spc."""
        spc = SPCFile(data_dir / "s_evenx.spc")
        assert spc.header["version"] == 0x4B
        assert spc.header["n_points"] == 1844
        assert len(spc) == 1
        assert spc.x.shape == (1844,)
        assert spc.y.shape == (1844,)
        assert pytest.approx(spc.x.min(), rel=1e-2) == 447.48
        assert pytest.approx(spc.x.max(), rel=1e-2) == 4002.28

    def test_m_evenz_properties(self, data_dir: Path) -> None:
        """Test known properties of m_evenz.spc."""
        spc = SPCFile(data_dir / "m_evenz.spc")
        assert spc.header["version"] == 0x4B
        assert len(spc) == 32
        assert spc.x.shape == (171,)
        assert spc.y.shape == (171, 32)
        assert pytest.approx(spc.x.min(), rel=1e-2) == 200.0
        assert pytest.approx(spc.x.max(), rel=1e-2) == 800.0


class TestUnitLabels:
    """Test human-readable unit label properties."""

    def test_unit_labels_single_spectrum(self, data_dir: Path) -> None:
        spc = SPCFile(data_dir / "s_evenx.spc")
        assert spc.x_unit == "Wavenumber (cm^-1)"
        assert spc.y_unit == "Absorbance (AU)"
        assert spc.z_unit == "Arbitrary"
        assert spc.w_unit == "Arbitrary"

    def test_unit_labels_explicit_x_file(self, data_dir: Path) -> None:
        spc = SPCFile(data_dir / "s_xy.spc")
        assert spc.x_unit == "Time (min)"
        assert spc.y_unit == "Arbitrary Intensity"

    def test_unit_labels_multifile(self, data_dir: Path) -> None:
        spc = SPCFile(data_dir / "m_evenz.spc")
        assert spc.x_unit == "Wavelength (nm)"
        assert spc.y_unit == "Absorbance (AU)"
        assert spc.z_unit == "Time (min)"

    def test_unit_labels_ftir(self, data_dir: Path) -> None:
        spc = SPCFile(data_dir / "ft-ir.spc")
        assert spc.x_unit == "Wavenumber (cm^-1)"
        assert spc.y_unit == "Transmittance"

    def test_unit_labels_mass_spec(self, data_dir: Path) -> None:
        spc = SPCFile(data_dir / "ms.spc")
        assert spc.x_unit == "Mass (m/z)"


class TestExperimentLabels:
    """Test human-readable experiment type label property."""

    @pytest.mark.parametrize(
        "filename",
        [
            "s_evenx.spc",
            "s_xy.spc",
            "m_evenz.spc",
            "ft-ir.spc",
            "ms.spc",
        ],
    )
    def test_experiment_label(self, data_dir: Path, filename: str) -> None:
        spc = SPCFile(data_dir / filename)
        assert spc.experiment == "General"


class TestDateField:
    """Test parsing of the packed date/time field."""

    @pytest.mark.parametrize("filename", ["s_evenx.spc", "ms.spc"])
    def test_date_none_for_missing_or_invalid(self, data_dir: Path, filename: str) -> None:
        spc = SPCFile(data_dir / filename)
        assert spc.date is None

    @pytest.mark.parametrize(
        ("filename", "expected"),
        [
            ("s_xy.spc", datetime(1986, 1, 9, 8, 47)),
            ("ft-ir.spc", datetime(1995, 4, 18, 9, 20)),
        ],
    )
    def test_date_known_values(self, data_dir: Path, filename: str, expected: datetime) -> None:
        spc = SPCFile(data_dir / filename)
        assert spc.date == expected


class TestHeaderFlagProperties:
    """Test convenience boolean properties for ftflgs bits."""

    @pytest.mark.parametrize(
        ("filename", "expected"),
        [
            (
                "s_evenx.spc",
                {
                    "y_16bit": False,
                    "chromatogram": False,
                    "multifile": False,
                    "random_z": False,
                    "ordered_z": False,
                    "custom_axis_labels": False,
                    "per_subfile_xy": False,
                    "explicit_x": False,
                },
            ),
            (
                "s_xy.spc",
                {
                    "y_16bit": False,
                    "chromatogram": False,
                    "multifile": False,
                    "random_z": False,
                    "ordered_z": False,
                    "custom_axis_labels": False,
                    "per_subfile_xy": False,
                    "explicit_x": True,
                },
            ),
            (
                "m_evenz.spc",
                {
                    "y_16bit": False,
                    "chromatogram": False,
                    "multifile": True,
                    "random_z": False,
                    "ordered_z": False,
                    "custom_axis_labels": False,
                    "per_subfile_xy": False,
                    "explicit_x": False,
                },
            ),
            (
                "ms.spc",
                {
                    "y_16bit": True,
                    "chromatogram": False,
                    "multifile": False,
                    "random_z": False,
                    "ordered_z": False,
                    "custom_axis_labels": True,
                    "per_subfile_xy": True,
                    "explicit_x": True,
                },
            ),
        ],
    )
    def test_flag_properties_match_known_files(self, data_dir: Path, filename: str, expected: dict[str, bool]) -> None:
        spc = SPCFile(data_dir / filename)

        assert spc.flags == int(spc.header["flags"])
        assert spc.y_16bit is expected["y_16bit"]
        assert spc.is_chromatogram is expected["chromatogram"]
        assert spc.is_multifile is expected["multifile"]
        assert spc.random_z is expected["random_z"]
        assert spc.ordered_z is expected["ordered_z"]
        assert spc.custom_axis_labels is expected["custom_axis_labels"]
        assert spc.per_subfile_xy is expected["per_subfile_xy"]
        assert spc.explicit_x is expected["explicit_x"]
