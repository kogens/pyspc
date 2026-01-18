import struct
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO

SPC_HEADER_SIZE = 512

SPC_HEADER_FIELDS: list[tuple[str, str]] = [
    ("flags", "B"),  # ftflgs: uint8
    ("version", "B"),  # fversn: uint8
    ("experiment_type", "B"),  # fexper: uint8
    ("exponent", "b"),  # fexp: int8
    ("n_points", "I"),  # fnpts: uint32
    ("first_x", "d"),  # ffirst: double
    ("last_x", "d"),  # flast: double
    ("n_subfiles", "I"),  # fnsub: uint32
    ("x_unit_code", "B"),  # fxtype: uint8
    ("y_unit_code", "B"),  # fytype: uint8
    ("z_unit_code", "B"),  # fztype: uint8
    ("posting_disposition", "B"),  # fpost: uint8
    ("date_time_raw", "I"),  # fdate: uint32
    ("resolution_str", "9s"),  # fres: char[9]
    ("source_str", "9s"),  # fsource: char[9]
    ("peak_point_index", "H"),  # fpeakpt: uint16
    ("spare_floats", "8f"),  # fspare: float[8]
    ("comment", "130s"),  # fcmnt: char[130]
    ("axis_label_text", "30s"),  # fcatxt: char[30]
    ("log_offset", "I"),  # flogoff: uint32
    ("modification_flags", "I"),  # fmods: uint32
    ("processing_code", "B"),  # fprocs: uint8
    ("calibration_level_raw", "B"),  # flevel: uint8
    ("sample_injection_number", "H"),  # fsampin: uint16
    ("data_multiplier", "f"),  # ffactor: float
    ("method_text", "48s"),  # fmethod: char[48]
    ("z_increment", "f"),  # fzinc: float
    ("w_planes", "I"),  # fwplanes: uint32
    ("w_increment", "f"),  # fwinc: float
    ("w_unit_code", "B"),  # fwtype: uint8
    ("reserved", "187s"),  # freserv: char[187]
]


def read_header(f: BinaryIO) -> dict[str, object]:
    # Ensure we are at the start of the file
    f.seek(0)
    buffer = f.read(SPC_HEADER_SIZE)
    if len(buffer) < SPC_HEADER_SIZE:
        raise ValueError(
            f"File too small to be a valid SPC file got {len(buffer)} bytes, expected {SPC_HEADER_SIZE} bytes"
        )
    header: dict[str, object] = {}
    offset = 0

    for field_name, fmt in SPC_HEADER_FIELDS:
        size = struct.calcsize(fmt)
        values = struct.unpack_from(fmt, buffer, offset)
        value = values[0] if len(values) == 1 else values

        # Interpret bytes as ASCII strings, stripping null bytes and whitespace
        if isinstance(value, bytes):
            value = value.split(b"\x00")[0].decode("latin-1", errors="replace").strip()

        header[field_name] = value
        offset += size

    return header
