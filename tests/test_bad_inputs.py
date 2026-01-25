from __future__ import annotations

from pathlib import Path

import pytest

from spcfile import SPCFile
from tests.fixtures_synthetic import choose_data_dir


@pytest.fixture(scope="session")
def data_dir(tmp_path_factory) -> Path:
    real = Path(__file__).parent / "data"
    return choose_data_dir(tmp_path_factory=tmp_path_factory, real_data_dir=real)


def test_truncated_header_raises_clear_value_error(tmp_path: Path, data_dir: Path) -> None:
    src = (data_dir / "s_evenx.spc").read_bytes()
    bad = tmp_path / "truncated_header.spc"
    bad.write_bytes(src[:100])

    with pytest.raises(ValueError, match=r"File too small: got 100 bytes, expected 512"):
        SPCFile(bad)


def test_truncated_subheader_raises_clear_value_error(tmp_path: Path, data_dir: Path) -> None:
    src = (data_dir / "s_evenx.spc").read_bytes()
    bad = tmp_path / "truncated_subheader.spc"
    # Keep a full header, then only part of the first SUBHDR.
    bad.write_bytes(src[: 512 + 10])

    with pytest.raises(ValueError, match=r"Could not read subheader: got 10 bytes, expected 32"):
        SPCFile(bad)


def test_truncated_explicit_x_payload_raises_clear_value_error(tmp_path: Path, data_dir: Path) -> None:
    src = (data_dir / "s_xy.spc").read_bytes()
    bad = tmp_path / "truncated_x.spc"

    # Truncate inside the global X array (this file uses TXVALS).
    bad.write_bytes(src[: 512 + 20])

    with pytest.raises(ValueError, match=r"Could not read explicit global X axis"):
        SPCFile(bad)


def test_truncated_y_payload_raises_clear_value_error(tmp_path: Path, data_dir: Path) -> None:
    src = (data_dir / "s_evenx.spc").read_bytes()
    bad = tmp_path / "truncated_y.spc"

    # Keep header + full SUBHDR, then truncate Y data.
    bad.write_bytes(src[: 512 + 32 + 100])

    with pytest.raises(ValueError, match=r"File ended while reading Y data"):
        SPCFile(bad)
