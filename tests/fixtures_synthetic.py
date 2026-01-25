from __future__ import annotations

import os
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

# This module builds small, deterministic SPC files for tests.
#
# Design goals:
# - Match the *current* reader implementation in spcfile/spcfile.py
# - Use the same filenames as tests/data so existing tests need minimal changes.
# - Keep layout stable: 0x4B little-endian new format only (except old-format stub).

SPC_HEADER_SIZE = 512
SPC_SUBHEADER_SIZE = 32

# Flag bits (duplicated from spcfile.spcfile to keep tests decoupled)
FLAG_Y_16BIT = 0x01
FLAG_IS_MULTIFILE = 0x04
FLAG_CUSTOM_AXIS_LABELS = 0x20
FLAG_PER_SUBFILE_XY = 0x40
FLAG_EXPLICIT_X = 0x80


SPC_HEADER_FIELDS: list[tuple[str, str]] = [
    ("flags", "B"),
    ("version", "B"),
    ("experiment_type", "B"),
    ("exponent", "b"),
    ("n_points", "I"),
    ("first_x", "d"),
    ("last_x", "d"),
    ("n_subfiles", "I"),
    ("x_unit_code", "B"),
    ("y_unit_code", "B"),
    ("z_unit_code", "B"),
    ("posting_disposition", "B"),
    ("date_int", "I"),
    ("resolution_str", "9s"),
    ("source_str", "9s"),
    ("peak_point_index", "H"),
    ("spare_floats", "8f"),
    ("comment", "130s"),
    ("axis_label_text", "30s"),
    ("log_offset", "I"),
    ("modification_flags", "I"),
    ("processing_code", "B"),
    ("calibration_level_raw", "B"),
    ("sample_injection_number", "H"),
    ("data_multiplier", "f"),
    ("method_text", "48s"),
    ("z_increment", "f"),
    ("w_planes", "I"),
    ("w_increment", "f"),
    ("w_unit_code", "B"),
    ("reserved", "187s"),
]

SPC_SUBHEADER_FIELDS: list[tuple[str, str]] = [
    ("flags", "B"),
    ("exponent", "b"),
    ("subfile_index", "H"),
    ("z_value", "f"),
    ("z_next", "f"),
    ("noise", "f"),
    ("n_points", "I"),
    ("n_scans", "I"),
    ("w_value", "f"),
    ("reserved", "4s"),
]


def _cstr(text: str, length: int) -> bytes:
    raw = text.encode("latin-1", errors="replace")
    raw = raw[:length]
    return raw + b"\x00" * (length - len(raw))


def _pack_date_int(year: int, month: int, day: int, hour: int, minute: int) -> int:
    # Packed: YYYY(12) MM(4) DD(5) HH(5) MM(6)
    return (year << 20) | (month << 16) | (day << 11) | (hour << 6) | minute


def _pack_header_le(fields: dict[str, object]) -> bytes:
    buf = bytearray()
    for name, fmt in SPC_HEADER_FIELDS:
        value = fields[name]
        if fmt.endswith("s"):
            size = int(fmt[:-1])
            assert isinstance(value, (bytes, bytearray))
            if len(value) != size:
                raise ValueError(f"Field {name} must be {size} bytes, got {len(value)}")
            buf += value
        elif fmt.endswith("f") and fmt != "f":
            # e.g. 8f
            assert isinstance(value, Iterable)
            buf += struct.pack("<" + fmt, *value)
        else:
            buf += struct.pack("<" + fmt, value)

    if len(buf) != SPC_HEADER_SIZE:
        raise ValueError(f"Header size mismatch: got {len(buf)} bytes, expected {SPC_HEADER_SIZE}")
    return bytes(buf)


def _pack_subheader_le(fields: dict[str, object]) -> bytes:
    buf = bytearray()
    for name, fmt in SPC_SUBHEADER_FIELDS:
        value = fields[name]
        if fmt.endswith("s"):
            size = int(fmt[:-1])
            assert isinstance(value, (bytes, bytearray))
            if len(value) != size:
                raise ValueError(f"Subheader field {name} must be {size} bytes, got {len(value)}")
            buf += value
        else:
            buf += struct.pack("<" + fmt, value)

    if len(buf) != SPC_SUBHEADER_SIZE:
        raise ValueError(f"Subheader size mismatch: got {len(buf)} bytes, expected {SPC_SUBHEADER_SIZE}")
    return bytes(buf)


def _encode_y_fixed_int32(y_float: np.ndarray, exponent: int) -> bytes:
    # Reader decodes: y = int32 * (2**exponent / 2**32)
    # So int32 = y * (2**32 / 2**exponent)
    scale = (2.0**32) / (2.0**exponent)
    y_raw = np.clip(np.round(y_float * scale), np.iinfo(np.int32).min, np.iinfo(np.int32).max).astype("<i4")
    return y_raw.tobytes(order="C")


def _encode_y_fixed_int16(y_float: np.ndarray, exponent: int) -> bytes:
    # Reader decodes: y = int16 * (2**exponent / 2**16)
    scale = (2.0**16) / (2.0**exponent)
    y_raw = np.clip(np.round(y_float * scale), np.iinfo(np.int16).min, np.iinfo(np.int16).max).astype("<i2")
    return y_raw.tobytes(order="C")


def _default_header(*, flags: int, n_points: int, n_subfiles: int) -> dict[str, object]:
    return {
        "flags": flags,
        "version": 0x4B,
        "experiment_type": 0,
        "exponent": 0,
        "n_points": int(n_points),
        "first_x": 0.0,
        "last_x": 1.0,
        "n_subfiles": int(n_subfiles),
        "x_unit_code": 0,
        "y_unit_code": 0,
        "z_unit_code": 0,
        "posting_disposition": 0,
        "date_int": 0,
        "resolution_str": _cstr("", 9),
        "source_str": _cstr("", 9),
        "peak_point_index": 0,
        "spare_floats": (0.0,) * 8,
        "comment": _cstr("", 130),
        "axis_label_text": _cstr("", 30),
        "log_offset": 0,
        "modification_flags": 0,
        "processing_code": 0,
        "calibration_level_raw": 0,
        "sample_injection_number": 0,
        "data_multiplier": 0.0,
        "method_text": _cstr("", 48),
        "z_increment": 0.0,
        "w_planes": 0,
        "w_increment": 0.0,
        "w_unit_code": 0,
        "reserved": b"\x00" * 187,
    }


def make_shared_x_single_evenx() -> bytes:
    n_points = 1844
    header = _default_header(flags=0, n_points=n_points, n_subfiles=1)
    header["first_x"] = 447.48
    header["last_x"] = 4002.28
    header["x_unit_code"] = 1  # Wavenumber
    header["y_unit_code"] = 2  # Absorbance

    rng = np.random.default_rng(123)
    y = rng.normal(loc=0.0, scale=0.2, size=n_points).astype(np.float64)

    subhdr = {
        "flags": 0,
        "exponent": 0,
        "subfile_index": 0,
        "z_value": 0.0,
        "z_next": 0.0,
        "noise": 0.0,
        "n_points": 0,
        "n_scans": 0,
        "w_value": 0.0,
        "reserved": b"\x00" * 4,
    }

    return b"".join([
        _pack_header_le(header),
        _pack_subheader_le(subhdr),
        _encode_y_fixed_int32(y, exponent=int(header["exponent"])),
    ])


def make_shared_x_single_explicit_x() -> bytes:
    n_points = 256
    flags = FLAG_EXPLICIT_X

    header = _default_header(flags=flags, n_points=n_points, n_subfiles=1)
    header["first_x"] = 0.0
    header["last_x"] = 10.0
    header["x_unit_code"] = 5  # Time (min)
    header["y_unit_code"] = 0  # Arbitrary Intensity
    header["date_int"] = _pack_date_int(1986, 1, 9, 8, 47)

    x_lin = np.linspace(header["first_x"], header["last_x"], n_points, dtype=np.float64)
    # Make it monotonic but not equal to linspace at rtol=1e-6.
    x = (x_lin + (np.arange(n_points, dtype=np.float64) / n_points) * 1e-3).astype("<f4")

    rng = np.random.default_rng(456)
    y = rng.normal(loc=0.0, scale=0.3, size=n_points).astype(np.float64)

    subhdr = {
        "flags": 0,
        "exponent": 0,
        "subfile_index": 0,
        "z_value": 0.0,
        "z_next": 0.0,
        "noise": 0.0,
        "n_points": 0,
        "n_scans": 0,
        "w_value": 0.0,
        "reserved": b"\x00" * 4,
    }

    return b"".join([
        _pack_header_le(header),
        x.tobytes(order="C"),
        _pack_subheader_le(subhdr),
        _encode_y_fixed_int32(y, exponent=int(header["exponent"])),
    ])


def make_shared_x_multifile_evenz() -> bytes:
    n_points = 171
    n_subfiles = 32
    flags = FLAG_IS_MULTIFILE

    header = _default_header(flags=flags, n_points=n_points, n_subfiles=n_subfiles)
    header["first_x"] = 200.0
    header["last_x"] = 800.0
    header["x_unit_code"] = 3  # Wavelength (nm)
    header["y_unit_code"] = 2  # Absorbance
    header["z_unit_code"] = 5  # Time (min)

    rng = np.random.default_rng(789)

    parts: list[bytes] = [_pack_header_le(header)]
    for i in range(n_subfiles):
        y = rng.normal(loc=0.0, scale=0.15, size=n_points).astype(np.float64)
        subhdr = {
            "flags": 0,
            "exponent": 0,
            "subfile_index": i,
            "z_value": float(i),
            "z_next": float(i + 1),
            "noise": 0.0,
            "n_points": 0,
            "n_scans": 0,
            "w_value": 0.0,
            "reserved": b"\x00" * 4,
        }
        parts.append(_pack_subheader_le(subhdr))
        parts.append(_encode_y_fixed_int32(y, exponent=int(header["exponent"])))

    return b"".join(parts)


def make_shared_x_4d_map() -> bytes:
    # Tests expect w_planes=11 and len=121
    n_points = 64
    n_subfiles = 121
    w_planes = 11
    flags = FLAG_IS_MULTIFILE

    header = _default_header(flags=flags, n_points=n_points, n_subfiles=n_subfiles)
    header["w_planes"] = w_planes
    header["x_unit_code"] = 0
    header["y_unit_code"] = 0

    rng = np.random.default_rng(2468)

    # 121 subfiles = 11 planes * 11 z positions
    z_per_plane = n_subfiles // w_planes

    parts: list[bytes] = [_pack_header_le(header)]
    for i in range(n_subfiles):
        plane = i // z_per_plane
        y = rng.normal(loc=0.0, scale=0.1, size=n_points).astype(np.float64)
        subhdr = {
            "flags": 0,
            "exponent": 0,
            "subfile_index": i,
            "z_value": float(i % z_per_plane),
            "z_next": float((i % z_per_plane) + 1),
            "noise": 0.0,
            "n_points": 0,
            "n_scans": 0,
            "w_value": float(plane),
            "reserved": b"\x00" * 4,
        }
        parts.append(_pack_subheader_le(subhdr))
        parts.append(_encode_y_fixed_int32(y, exponent=int(header["exponent"])))

    return b"".join(parts)


def make_log_text_block(text: str) -> bytes:
    # Reader expects 64-byte LOGSTC header; it unpacks first 20 bytes as 5 uint32.
    text_bytes = text.encode("latin-1", errors="replace")
    txt_offset = 64
    logsize = 64 + len(text_bytes) + 1  # include NUL terminator

    logstc = bytearray(64)
    struct.pack_into("<IIIII", logstc, 0, logsize, logsize, txt_offset, 0, 0)
    # Remaining 44 bytes already zero.

    return bytes(logstc) + text_bytes + b"\x00"


def make_shared_x_single_with_log() -> bytes:
    # Modelled after ft-ir.spc tests.
    n_points = 256
    header = _default_header(flags=0, n_points=n_points, n_subfiles=1)
    header["first_x"] = 447.48
    header["last_x"] = 4002.28
    header["x_unit_code"] = 1  # Wavenumber
    header["y_unit_code"] = 128  # Transmittance
    header["date_int"] = _pack_date_int(1995, 4, 18, 9, 20)

    rng = np.random.default_rng(1357)
    y = rng.normal(loc=0.0, scale=0.12, size=n_points).astype(np.float64)

    subhdr = {
        "flags": 0,
        "exponent": 0,
        "subfile_index": 0,
        "z_value": 0.0,
        "z_next": 0.0,
        "noise": 0.0,
        "n_points": 0,
        "n_scans": 0,
        "w_value": 0.0,
        "reserved": b"\x00" * 4,
    }

    # Build payload without log first so we can set log_offset.
    payload = b"".join([
        _pack_header_le(header),
        _pack_subheader_le(subhdr),
        _encode_y_fixed_int32(y, exponent=int(header["exponent"])),
    ])

    # Create deterministic log text with exact length expected by tests.
    # Tests assert startswith("MODEL") and len == 376.
    base = "MODEL=synthetic\r\nSCANS=1\r\n"
    if len(base) > 376:
        raise ValueError("Base log text unexpectedly long")
    log_text = base + ("X" * (376 - len(base)))

    log_block = make_log_text_block(log_text)

    # Patch log_offset in header.
    log_offset = len(payload)
    header["log_offset"] = log_offset

    payload_with_log = b"".join([
        _pack_header_le(header),
        _pack_subheader_le(subhdr),
        _encode_y_fixed_int32(y, exponent=int(header["exponent"])),
        log_block,
    ])

    return payload_with_log


def make_xyxy_single_ms() -> bytes:
    # Single-subfile TXYXYS set: reader will still expose spc.x/spc.y.
    n_points = 50
    flags = FLAG_Y_16BIT | FLAG_CUSTOM_AXIS_LABELS | FLAG_PER_SUBFILE_XY | FLAG_EXPLICIT_X

    header = _default_header(flags=flags, n_points=0, n_subfiles=1)
    header["x_unit_code"] = 9  # Mass (m/z)
    header["y_unit_code"] = 0
    header["axis_label_text"] = _cstr("X\x00Y\x00", 30)

    x = np.linspace(10.0, 100.0, n_points, dtype=np.float64).astype("<f4")
    rng = np.random.default_rng(4242)
    y = rng.normal(loc=0.0, scale=0.2, size=n_points).astype(np.float64)

    subhdr = {
        "flags": 0,
        "exponent": 0,
        "subfile_index": 0,
        "z_value": 0.0,
        "z_next": 0.0,
        "noise": 0.0,
        "n_points": n_points,
        "n_scans": 0,
        "w_value": 0.0,
        "reserved": b"\x00" * 4,
    }

    return b"".join([
        _pack_header_le(header),
        _pack_subheader_le(subhdr),
        x.tobytes(order="C"),
        _encode_y_fixed_int16(y, exponent=int(header["exponent"])),
    ])


def make_xyxy_multifile_m_xyxy() -> bytes:
    flags = FLAG_IS_MULTIFILE | FLAG_PER_SUBFILE_XY | FLAG_EXPLICIT_X

    n_subfiles = 3
    header = _default_header(flags=flags, n_points=0, n_subfiles=n_subfiles)
    header["x_unit_code"] = 0
    header["y_unit_code"] = 0

    rng = np.random.default_rng(999)

    parts: list[bytes] = [_pack_header_le(header)]
    for i, n_points in enumerate((10, 12, 8)):
        x = np.linspace(1.0, float(n_points), n_points, dtype=np.float64).astype("<f4")
        y = rng.normal(loc=0.0, scale=0.1, size=n_points).astype(np.float64)

        subhdr = {
            "flags": 0,
            "exponent": 0,
            "subfile_index": i,
            "z_value": float(i),
            "z_next": float(i),
            "noise": 0.0,
            "n_points": n_points,
            "n_scans": 0,
            "w_value": 0.0,
            "reserved": b"\x00" * 4,
        }

        parts.append(_pack_subheader_le(subhdr))
        parts.append(x.tobytes(order="C"))
        parts.append(_encode_y_fixed_int32(y, exponent=int(header["exponent"])))

    return b"".join(parts)


def make_old_format_stub() -> bytes:
    # The reader rejects version 0x4D before parsing. We only need a 512-byte buffer
    # with byte[1] = 0x4D.
    buf = bytearray(b"\x00" * 512)
    buf[1] = 0x4D
    return bytes(buf)


@dataclass(frozen=True)
class SyntheticDataset:
    """A synthetic dataset that mimics tests/data filenames."""

    files: dict[str, bytes]


def build_synthetic_dataset() -> SyntheticDataset:
    return SyntheticDataset(
        files={
            "s_evenx.spc": make_shared_x_single_evenx(),
            "s_xy.spc": make_shared_x_single_explicit_x(),
            "m_evenz.spc": make_shared_x_multifile_evenz(),
            "nir.spc": make_shared_x_multifile_evenz(),
            "4d_map.spc": make_shared_x_4d_map(),
            "ft-ir.spc": make_shared_x_single_with_log(),
            "raman.spc": make_shared_x_single_with_log(),
            "nmr_fid.spc": make_shared_x_single_with_log(),
            "nmr_spc.spc": make_shared_x_single_with_log(),
            "ms.spc": make_xyxy_single_ms(),
            "m_xyxy.spc": make_xyxy_multifile_m_xyxy(),
            "m_ordz.spc": make_old_format_stub(),
        }
    )


def write_synthetic_dataset(dir_path: Path) -> Path:
    dir_path.mkdir(parents=True, exist_ok=True)
    dataset = build_synthetic_dataset()
    for name, payload in dataset.files.items():
        (dir_path / name).write_bytes(payload)

    # Keep the directory looking like tests/data.
    (dir_path / "testdata.md").write_text(
        "Synthetic dataset generated at test time.\n",
        encoding="utf-8",
    )

    return dir_path


def choose_data_dir(*, tmp_path_factory, real_data_dir: Path) -> Path:
    """Return a directory containing .spc files for tests.

    Preference order:
    1) Use real_data_dir if it looks populated.
    2) Otherwise generate a synthetic dataset into a temp directory.

    You can force synthetic fixtures by setting SPC_USE_SYNTHETIC_DATA=1.
    """

    force = os.environ.get("SPC_USE_SYNTHETIC_DATA", "").strip() == "1"

    # "Looks populated" means at least a couple of key files exist.
    has_real = (real_data_dir / "s_evenx.spc").is_file() and (real_data_dir / "m_evenz.spc").is_file()

    if has_real and not force:
        return real_data_dir

    synth_dir = tmp_path_factory.mktemp("spc_synthetic_data")
    return write_synthetic_dataset(Path(synth_dir))
