from pathlib import Path
import numpy as np
import pytest

from pyspc import SPCFile, SPCSubfile
from pyspc.spcfile import FLAG_EXPLICIT_X


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return Path(__file__).parent / "data"


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

    def test_4d_file_not_supported(self, data_dir: Path) -> None:
        """Raise NotImplementedError for 4D W-plane files until supported."""
        with pytest.raises(NotImplementedError):
            SPCFile(data_dir / "4d_map.spc")

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
        assert len(text)  == 376  # Exact log length for this file

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
