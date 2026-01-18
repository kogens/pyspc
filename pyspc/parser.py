"""Low-level binary parsing for SPC files."""

import struct
from typing import BinaryIO

import numpy as np

SPC_HEADER_SIZE = 512
SPC_SUBHEADER_SIZE = 32

# Main header flag bits (ftflgs)
FLAG_Y_16BIT = 0x01  # TSPREC: Y data is 16-bit if set, 32-bit if clear
FLAG_CHROMATOGRAM = 0x02  # TCGRAM: Enables chromatogram / fexper interpretation
FLAG_MULTIFILE = 0x04  # TMULTI: Multifile; more than one subfile present
FLAG_RANDOM_Z = 0x08  # TRANDM: Multifile with arbitrary (unordered) Z values
FLAG_ORDERED_Z = 0x10  # TORDRD: Multifile with ordered but uneven Z values
FLAG_CUSTOM_AXIS_LABELS = 0x20  # TALABS: Use fcatxt axis labels instead of type defaults
FLAG_PER_SUBFILE_XY = 0x40  # TXYXYS: Per-subfile X arrays and lengths (TXYXYS mode)
FLAG_EXPLICIT_X = 0x80  # TXVALS: X values stored explicitly as float array(s)


SPC_HEADER_FIELDS: list[tuple[str, str]] = [
    ("flags", "B"),  # ftflgs: uint8, file flags
    ("version", "B"),  # fversn: uint8, SPC format version
    ("experiment_type", "B"),  # fexper: uint8, experiment type
    ("exponent", "b"),  # fexp: int8, data exponent
    ("n_points", "I"),  # fnpts: uint32, number of points per spectrum
    ("first_x", "d"),  # ffirst: double, first x value
    ("last_x", "d"),  # flast: double, last x value
    ("n_subfiles", "I"),  # fnsub: uint32, number of subfilesTSPREC)
    ("x_unit_code", "B"),  # fxtype: uint8, x axis unit code
    ("y_unit_code", "B"),  # fytype: uint8, y axis unit code
    ("z_unit_code", "B"),  # fztype: uint8, z axis unit code
    ("posting_disposition", "B"),  # fpost: uint8, posting disposition ()
    ("date_time_raw", "I"),  # fdate: uint32, datetime as packed integer
    (
        "resolution_str",
        "9s",
    ),  # fres: char[9], resolution string (e.g. "4 cm-1" per data point)
    ("source_str", "9s"),  # fsource: char[9], instrument or method
    ("peak_point_index", "H"),  # fpeakpt: uint16, index of peak point
    ("spare_floats", "8f"),  # fspare: float[8], unused, reserved for future use
    ("comment", "130s"),  # fcmnt: char[130], comment string
    (
        "axis_label_text",
        "30s",
    ),  # fcatxt: char[30], x axis label text (used if TALABS flag is set)
    (
        "log_offset",
        "I",
    ),  # flogoff: uint32, file offset (in bytes) to start of log data (0 = no log)
    (
        "modification_flags",
        "I",
    ),  # fmods: uint32, bitfield indicating file modifications (used by GRAMS)
    (
        "processing_code",
        "B",
    ),  # fprocs: uint8, code indicating type of processing done to data
    (
        "calibration_level_raw",
        "B",
    ),  # flevel: uint8, calibration status of data (1=not calibration, >1 = calibration data)
    (
        "sample_injection_number",
        "H",
    ),  # fsampin: uint16, sample injection number for sub-methods
    ("data_multiplier", "f"),  # ffactor: float, data multiplier factor for y values
    ("method_text", "48s"),  # fmethod: char[48], method or technique text
    ("z_increment", "f"),  # fzinc: float, z axis increment
    ("w_planes", "I"),  # fwplanes: uint32, number of w planes (3rd dimension)
    ("w_increment", "f"),  # fwinc: float, w axis increment
    ("w_unit_code", "B"),  # fwtype: uint8, w axis unit code
    ("reserved", "187s"),  # freserv: char[187], unused
]

SPC_SUBHEADER_FIELDS: list[tuple[str, str]] = [
    ("flags", "B"),  # subflgs: uint8, flags for this subfile
    ("exponent", "b"),  # subexp: int8, exponent for this subfile
    ("subfile_index", "H"),  # subindx: uint16, index of this subfile
    ("z_value", "f"),  # subtime: float, z value for this subfile
    ("z_next", "f"),  # subnext: float, z value of next subfile
    ("noise", "f"),  # subnois: float, noise level for interferogram
    ("n_points", "I"),  # subnpts: uint32, number of datapoints in this subfile
    ("n_scans", "I"),  # subscan: uint32, number of scans averaged
    ("w_value", "f"),  # subwlevel: float, w value for this subfile
    ("reserved", "4s"),  # subresv: char[4], reserved (unused)
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


def read_subheader(f: BinaryIO, offset: int) -> dict[str, object]:
    """Read an SPC subheader from the given file at the specified offset."""
    f.seek(offset)
    buffer = f.read(SPC_SUBHEADER_SIZE)
    if len(buffer) < SPC_SUBHEADER_SIZE:
        raise ValueError(
            f"Could not read full subheader at offset {offset}, got {len(buffer)} bytes, expected {SPC_SUBHEADER_SIZE} bytes"
        )
    subheader: dict[str, object] = {}
    sub_offset = 0

    for field_name, fmt in SPC_SUBHEADER_FIELDS:
        size = struct.calcsize(fmt)
        values = struct.unpack_from(fmt, buffer, sub_offset)
        value = values[0] if len(values) == 1 else values  #

        subheader[field_name] = value
        sub_offset += size
    return subheader


def read_x_array(
    f: BinaryIO,
    n_points: int,
    flags: int,
    has_global_x: bool,
    global_x: np.ndarray | None,
    first_x: float,
    last_x: float,
) -> np.ndarray:
    """Read X coordinate array for a subfile.

    Args:
        f: Open binary file at the correct position
        n_points: Number of points to read
        flags: Main header flags
        has_global_x: Whether a global X array was already read
        global_x: The global X array if present
        first_x: First X value from main header (for evenly spaced)
        last_x: Last X value from main header (for evenly spaced)

    Returns:
        numpy array of X coordinates
    """
    if flags & FLAG_PER_SUBFILE_XY:
        # Read per-subfile X array
        return np.frombuffer(f.read(n_points * 4), dtype="<f4")
    elif has_global_x:
        # Use global X array
        return global_x
    else:
        # Evenly spaced X - compute from header
        return np.linspace(first_x, last_x, n_points)


def read_y_array(
    f: BinaryIO,
    n_points: int,
    y_exponent: int,
    is_16bit: bool,
) -> np.ndarray:
    """Read and decode Y value array for a subfile.

    Args:
        f: Open binary file at the correct position
        n_points: Number of points to read
        y_exponent: Exponent for fixed-point scaling (or -128 for float)
        is_16bit: Whether Y values are 16-bit (else 32-bit)

    Returns:
        numpy array of Y values as float64
    """
    # Read raw Y data
    if is_16bit:
        y_raw = np.frombuffer(f.read(n_points * 2), dtype="<i2")
    else:
        y_raw = np.frombuffer(f.read(n_points * 4), dtype="<i4")

    # Convert Y to float
    if y_exponent == -128:  # 0x80 = floating point Y
        # Re-interpret raw bytes as float32
        if is_16bit:
            raise ValueError("Cannot have 16-bit Y with float exponent")
        y_data = y_raw.view("<f4").astype(np.float64)
    else:
        # Fixed-point conversion
        bit_width = 16 if is_16bit else 32
        scale = 2.0**y_exponent / (2.0**bit_width)
        y_data = y_raw.astype(np.float64) * scale

    return y_data


def read_subfile(
    f: BinaryIO,
    main_header: dict[str, object],
    global_x: np.ndarray | None,
    has_global_x: bool,
    is_16bit_y: bool,
) -> dict[str, object]:
    """Read one complete subfile (header + X + Y data).

    Args:
        f: Open binary file positioned at the start of a subfile header
        main_header: Parsed main header
        global_x: Global X array if present, else None
        has_global_x: Whether global X array exists
        is_16bit_y: Whether Y values are 16-bit

    Returns:
        Dict with 'header', 'x', and 'y' keys
    """
    flags = main_header["flags"]
    n_points = main_header["n_points"]
    exponent = main_header["exponent"]

    # Read and parse subheader
    subhdr_buffer = f.read(SPC_SUBHEADER_SIZE)
    if len(subhdr_buffer) < SPC_SUBHEADER_SIZE:
        raise ValueError(
            f"Could not read full subheader, got {len(subhdr_buffer)} bytes"
        )

    subheader: dict[str, object] = {}
    offset = 0
    for field_name, fmt in SPC_SUBHEADER_FIELDS:
        size = struct.calcsize(fmt)
        values = struct.unpack_from(fmt, subhdr_buffer, offset)
        value = values[0] if len(values) == 1 else values
        subheader[field_name] = value
        offset += size

    # Determine number of points for this subfile
    if flags & FLAG_PER_SUBFILE_XY:
        points_in_subfile = subheader["n_points"]
    else:
        points_in_subfile = n_points

    # Read X array
    x_data = read_x_array(
        f,
        points_in_subfile,
        flags,
        has_global_x,
        global_x,
        main_header["first_x"],
        main_header["last_x"],
    )

    # Determine Y exponent for this subfile
    sub_exponent = subheader["exponent"]
    if flags & FLAG_MULTIFILE:
        y_exponent = sub_exponent
    else:
        y_exponent = exponent

    # Read Y array
    y_data = read_y_array(f, points_in_subfile, y_exponent, is_16bit_y)

    return {"header": subheader, "x": x_data, "y": y_data}


def read_all_subfiles(
    f: BinaryIO, header: dict[str, object]
) -> list[dict[str, object]]:
    """Read all complete subfiles (header + X/Y data) from an SPC file.

    Args:
        f: Open binary file positioned at start (will seek as needed)
        header: Parsed main header dict from read_header()

    Returns:
        List of subfile dicts, each containing:
        - 'header': subheader dict
        - 'x': numpy array of X values
        - 'y': numpy array of Y values

    Reads each subfile completely in a single sequential pass through the file.
    """
    flags = header["flags"]
    n_points = header["n_points"]
    n_subfiles = header["n_subfiles"]

    # Start after main header
    f.seek(SPC_HEADER_SIZE)

    # Read global X array if present (FLAG_EXPLICIT_X set, FLAG_PER_SUBFILE_XY clear)
    global_x = None
    has_global_x = (flags & FLAG_EXPLICIT_X) and not (flags & FLAG_PER_SUBFILE_XY)
    if has_global_x:
        global_x = np.frombuffer(f.read(n_points * 4), dtype="<f4")

    # Determine Y storage format
    is_16bit_y = bool(flags & FLAG_Y_16BIT)

    # Read each complete subfile
    subfiles: list[dict[str, object]] = []
    for i in range(n_subfiles):
        subfile = read_subfile(f, header, global_x, has_global_x, is_16bit_y)
        subfiles.append(subfile)

    return subfiles
