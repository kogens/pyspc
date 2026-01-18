from pathlib import Path
from pyspc import SPCFile, SPCSubfile
from pyspc.parser import read_header

import pytest


@pytest.fixture
def data_dir() -> Path:
    return Path(__file__).parent / "data"

# def test_s_evenx_basic(data_dir: Path) -> None:
#     # Test loading a basic SPC file; single spectrum with even X spacing
#     spc = SPCFile(data_dir / "s_evenx.spc")
#     print(spc.header)
#     assert spc.x.ndim == 1
#     assert spc.y.ndim == 1
#     assert spc.x.size == spc.y.size > 0


class TestSPCFile:
    def test_s_evenx_basic(self, data_dir: Path) -> None:
        # Test loading a basic SPC file; single spectrum with even X spacing
        spc = SPCFile(data_dir / "s_evenx.spc")
        assert spc.x.ndim == 1
        assert spc.y.ndim == 1
        assert spc.x.size == spc.y.size > 0

class TestSPCHeader:
    def test_header_values_s_evenx(self, data_dir: Path) -> None:
        with open(data_dir / "s_evenx.spc", "rb") as f:
            header = read_header(f)
        assert header["version"] == 0x4B
        assert header["experiment_type"] == 0
        assert header['n_points'] == 1844