"""GRAMS SPC file reader."""

from __future__ import annotations

import datetime
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

# Constants and lookup tables
SPC_HEADER_SIZE = 512
SPC_SUBHEADER_SIZE = 32

# Main header flag bits (ftflgs)
FLAG_Y_16BIT = 0x01  # Y data stored as 16‑bit integers if set, 32‑bit if clear (TSPREC)
FLAG_CHROMATOGRAM = 0x02  # Chromatogram / special fexper interpretation (TCGRAM)
FLAG_IS_MULTIFILE = 0x04  # File contains multiple subfiles/spectra (TMULTI)
FLAG_RANDOM_Z = 0x08  # Multifile with arbitrary, unordered Z values (TRANDM)
FLAG_ORDERED_Z = 0x10  # Multifile with ordered but unevenly spaced Z values (TORDRD)
FLAG_CUSTOM_AXIS_LABELS = 0x20  # Use custom axis label text from header instead of defaults (TALABS)
FLAG_PER_SUBFILE_XY = 0x40  # Each subfile has its own X grid and length (XYXY / TXYXYS)
FLAG_EXPLICIT_X = 0x80  # A shared X axis is stored explicitly as float array(s) (TXVALS)


SPC_EXPERIMENT_TYPES: dict[int, str] = {
    0: "General",
    1: "Gas chromatogram",
    2: "Liquid chromatogram",
    3: "FT-IR",
    4: "NIR",
    5: "UV-VIS",
    6: "X-ray diffraction",
    7: "Mass spectrum",
    8: "NMR",
    9: "Raman",
    10: "Fluorescence",
    11: "Atomic",
    12: "Chromatogram (general)",
    13: "Color",
    14: "Simulated",
}

X_UNIT_LABELS: dict[int, str] = {
    0: "Arbitrary",
    1: "Wavenumber (cm^-1)",
    2: "Wavelength (µm)",
    3: "Wavelength (nm)",
    4: "Time (s)",
    5: "Time (min)",
    6: "Frequency (Hz)",
    7: "Frequency (kHz)",
    8: "Frequency (MHz)",
    9: "Mass (m/z)",
    10: "Parts per million (ppm)",
    11: "Time (days)",
    12: "Time (years)",
    13: "Raman Shift (cm^-1)",
    14: "Energy (eV)",
    15: "XYZ Text",
    16: "Diode number",
    17: "Channel",
    18: "Angle (deg)",
    19: "Temperature (F)",
    20: "Temperature (C)",
    21: "Temperature (K)",
    22: "Data points",
    23: "Time (ms)",
    24: "Time (µs)",
    25: "Time (ns)",
    26: "Frequency (GHz)",
    27: "Distance (cm)",
    28: "Distance (m)",
    29: "Distance (mm)",
    30: "Time (hours)",
    255: "No Units",
}

Y_UNIT_LABELS: dict[int, str] = {
    0: "Arbitrary Intensity",
    1: "Interferogram",
    2: "Absorbance (AU)",
    3: "Kubelka-Munk",
    4: "Counts",
    5: "Voltage (V)",
    6: "Angle (deg)",
    7: "Current (mA)",
    8: "Distance (mm)",
    9: "Voltage (mV)",
    10: "log(1/R)",
    11: "Percent (%)",
    12: "Intensity",
    13: "Relative Intensity",
    14: "Energy",
    16: "Decibel (dB)",
    19: "Temperature (F)",
    20: "Temperature (C)",
    21: "Temperature (K)",
    22: "Index of refraction",
    23: "Extinction coeff.",
    24: "Real",
    25: "Imaginary",
    26: "Complex",
    128: "Transmittance",
    129: "Reflectance",
    130: "Valley",
    255: "No units",
}

# Define the main header and subheader fields and their struct formats for parsing
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


@dataclass
class SPCSubfile:
    """Single spectrum extracted from an SPC file."""

    x: np.ndarray
    y: np.ndarray
    z: float | None = None
    subheader: dict[str, object] | None = None

    def __repr__(self) -> str:
        return f"<SPCSubfile n_points={len(self.x)} z={self.z}>"

    def __getattr__(self, name) -> object:
        """Returns the attributes as keys from the subheader dictionary if it exists."""
        if self.subheader and name in self.subheader:
            return self.subheader[name]
        raise AttributeError(f"'SPCSubfile' object has no attribute '{name}'")


class SPCFile:
    """In-memory representation of a GRAMS SPC file.

    Primary data:
        x: 1D array of X coordinates (shared across all subfiles)
        y: 1D array for single spectrum, 2D [n_points, n_subfiles] for multifile
        header: Main file header dict
        subheaders: List of per-subfile header dicts

    Indexing:
        spc[k] returns SPCSubfile with x, y for k-th spectrum
    """

    def __init__(self, path: str | Path) -> None:
        self._path = Path(path)

        if not self._path.is_file():
            raise FileNotFoundError(self._path)

        with self._path.open("rb") as f:
            header_buf = self._read_header_bytes(f)

            version = int(header_buf[1])
            if version == 0x4B:
                self._struct_prefix = "<"  # little-endian new format
            elif version == 0x4C:
                self._struct_prefix = ">"  # big-endian new format
            else:
                raise ValueError(f"Unsupported SPC version: {version:02X}")

            self.header = self._parse_header(header_buf, self._struct_prefix)

            if int(self.header.get("w_planes", 0)) < 0:
                raise ValueError(f"Invalid w_planes: {self.header.get('w_planes')}")

            if self.per_subfile_xy:
                # XYXY mode: each subfile has its own X grid and point count
                self.subfiles, self.subheaders = self._read_all_subfiles_xyxys(f)

                # If there is only one subfile, expose its X/Y directly via .x/.y
                # for convenience (even though TXYXYS is set).
                if not self.is_multifile:
                    self._x = self.subfiles[0].x
                    self._y = self.subfiles[0].y
                else:
                    self._x = np.array([], dtype=np.float64)
                    self._y = np.array([], dtype=np.float64)
            else:
                # Read X axis
                self._x = self._read_x_axis(f)

                # Read all Y data and subheaders
                self._y, self.subheaders = self._read_all_subfiles(f)
                self.subfiles = self._build_subfiles_shared_x(self._x, self._y, self.subheaders)

            # Read log block
            self.log = self._read_log_block(f)

    @property
    def flags(self) -> int:
        """Raw main-header flags bitfield (ftflgs)."""
        return int(self.header["flags"])

    @property
    def z(self) -> np.ndarray:
        """Z axis values.

        - Multifile, shared-X: 1D array of z_value per subfile.
        - Single-spectrum (including single XYXY): scalar Z coordinate (or 0.0).
        """
        z_values = [subheader.get("z_value", 0.0) for subheader in self.subheaders]
        return np.array(z_values)

    @property
    def w_planes(self) -> int:
        """Number of W-planes (0 for non-4D files)."""
        return int(self.header.get("w_planes", 0))

    @property
    def w(self) -> np.ndarray:
        """W coordinate per subfile.

        Returns a 1D array of ``w_value`` from each SUBHDR, one per subfile.

        For 4D data, values will repeat in blocks corresponding to Z positions
        within each W-plane; ``w_planes`` still exposes the plane count.
        """
        w_values = [float(subheader.get("w_value", 0.0)) for subheader in self.subheaders]
        return np.array(w_values)

    @property
    def path(self) -> Path:
        """Filesystem path of the underlying SPC file."""
        return self._path

    @staticmethod
    def _map_code(code: object, labels: dict[int, str]) -> str:
        code_int = int(code)
        return labels.get(code_int, f"Unknown ({code_int})")

    @property
    def experiment(self) -> str:
        """Human-readable experiment type."""
        return self._map_code(self.header["experiment_type"], SPC_EXPERIMENT_TYPES)

    @property
    def x_unit(self) -> str:
        """Human-readable X axis unit label."""
        return self._map_code(self.header["x_unit_code"], X_UNIT_LABELS)

    @property
    def y_unit(self) -> str:
        """Human-readable Y axis unit label."""
        return self._map_code(self.header["y_unit_code"], Y_UNIT_LABELS)

    @property
    def z_unit(self) -> str:
        """Human-readable Z axis unit label."""
        return self._map_code(self.header["z_unit_code"], X_UNIT_LABELS)

    @property
    def w_unit(self) -> str:
        """Human-readable W axis unit label."""
        return self._map_code(self.header["w_unit_code"], X_UNIT_LABELS)

    @property
    def date(self) -> datetime.datetime | None:
        """Date and time string from the header. Interpreted as packed int; YYYY(20) MM(4) DD(5)"""
        dt_raw = self.header["date_int"]
        minute = dt_raw & 0x3F
        hour = (dt_raw >> 6) & 0x1F
        day = (dt_raw >> 11) & 0x1F
        month = (dt_raw >> 16) & 0x0F
        year = (dt_raw >> 20) & 0xFFF

        # Basic validation
        if not ((1 <= month <= 12) and (1 <= day <= 31) and (1900 < year < 2100)):
            return None

        return datetime.datetime(year, month, day, hour, minute)

    @property
    def has_shared_x(self) -> bool:
        """True if all subfiles share a common X axis.

        False for XYXY (TXYXYS) files where each subfile has its own X array.
        """
        return (not self.per_subfile_xy) or (not self.is_multifile)

    @property
    def y_16bit(self) -> bool:
        """True if Y values are stored as 16-bit words (TSPREC)."""
        return bool(self.flags & FLAG_Y_16BIT)

    @property
    def is_chromatogram(self) -> bool:
        """True if chromatogram semantics are enabled (TCGRAM)."""
        return bool(self.flags & FLAG_CHROMATOGRAM)

    @property
    def is_multifile(self) -> bool:
        """True if the file contains multiple subfiles (TMULTI)."""
        return bool(self.flags & FLAG_IS_MULTIFILE)

    @property
    def random_z(self) -> bool:
        """True if multifile Z values are arbitrary/unordered (TRANDM)."""
        return bool(self.flags & FLAG_RANDOM_Z)

    @property
    def ordered_z(self) -> bool:
        """True if multifile Z values are ordered but uneven (TORDRD)."""
        return bool(self.flags & FLAG_ORDERED_Z)

    @property
    def custom_axis_labels(self) -> bool:
        """True if fcatxt custom axis labels are used (TALABS)."""
        return bool(self.flags & FLAG_CUSTOM_AXIS_LABELS)

    @property
    def per_subfile_xy(self) -> bool:
        """True if each subfile has its own X grid and length (XYXY / TXYXYS)."""
        return bool(self.flags & FLAG_PER_SUBFILE_XY)

    @property
    def explicit_x(self) -> bool:
        """True if a shared X axis is stored explicitly as float array(s) (TXVALS)."""
        return bool(self.flags & FLAG_EXPLICIT_X)

    @property
    def x(self) -> np.ndarray:
        """X coordinates (only for files with shared X axis).

        Returns:
            1D array of X coordinates shared across all subfiles.

        Raises:
            ValueError: For XYXY (TXYXYS) files with per-subfile X arrays.
                        Use spc[i].x to access the X grid for each spectrum.
        """
        if not self.is_multifile:
            return self.subfiles[0].x

        if self.is_multifile and not self.has_shared_x:
            raise ValueError(
                "This SPC file uses per-subfile X arrays (XYXY / TXYXYS). Use spc[i].x to access X for each spectrum."
            )
        return self._x

    @property
    def y(self) -> np.ndarray:
        """Y values (only for files with shared X axis).

        Returns:
            1D array for single spectrum.
            2D array [n_points, n_subfiles] for multifile.

        Raises:
            ValueError: For XYXY (TXYXYS) files with varying lengths.
                        Use spc[i].y to access individual Y arrays.
        """
        if not self.is_multifile:
            return self.subfiles[0].y

        if not self.has_shared_x:
            raise ValueError(
                "This SPC file uses per-subfile XY arrays (XYXY / TXYXYS). Use spc[i].y to access Y for each spectrum."
            )
        return self._y

    def __len__(self) -> int:
        """Number of subfiles (spectra) in the file."""
        return len(self.subfiles)

    def __getitem__(self, index: int) -> SPCSubfile:
        """Get k-th spectrum as an SPCSubfile."""
        return self.subfiles[index]

    def __getattr__(self, name: str) -> object:
        """Returns the attributes as keys from the header dictionary if it exists."""
        if name in self.header:
            return self.header[name]
        raise AttributeError(f"'SPCFile' object has no attribute '{name}'")

    def __iter__(self) -> Iterable[SPCSubfile]:
        """Iterate over all subfiles."""
        for i in range(len(self)):
            yield self[i]

    def __repr__(self) -> str:
        flags = self.header["flags"]
        parts = ["multifile" if flags & FLAG_IS_MULTIFILE else "single"]
        parts.append("shared-x" if self.has_shared_x else "per-subfile-x")
        if flags & FLAG_EXPLICIT_X:
            parts.append("explicit-x")
        flags = "|".join(parts)

        return f"<SPCFile path={self._path!r} n_subfiles={len(self)} n_points={self.header['n_points']} flags={flags}>"

    def __str__(self) -> str:
        if self.has_shared_x:
            pts = f"{self.header['n_points']}"
        else:
            pts = "Variable (per‑subfile X/Y)"

        lines = [
            f"SPC File: {self._path}",
            f"Date: {self.date}",
            f"Subfiles: {len(self)}",
            f"Points per subfile: {pts}",
            f"Experiment type: {self.experiment}",
            f"Units: X='{self.x_unit}', Y='{self.y_unit}', Z='{self.z_unit}', W='{self.w_unit}'",
        ]
        return "\n".join(lines)

    def _read_header_bytes(self, f) -> bytes:
        """Read the raw 512-byte main header bytes."""
        f.seek(0)
        buffer = f.read(SPC_HEADER_SIZE)
        if len(buffer) < SPC_HEADER_SIZE:
            raise ValueError(f"File too small: got {len(buffer)} bytes, expected {SPC_HEADER_SIZE}")
        return buffer

    @staticmethod
    def _parse_header(buffer: bytes, struct_prefix: str) -> dict[str, object]:
        """Parse the 512-byte main header.

        Args:
            buffer: Raw 512-byte SPCHDR.
            struct_prefix: '<' for little-endian (0x4B) or '>' for big-endian (0x4C).
        """
        if len(buffer) < SPC_HEADER_SIZE:
            raise ValueError(f"File too small: got {len(buffer)} bytes, expected {SPC_HEADER_SIZE}")

        header: dict[str, object] = {}
        offset = 0

        for field_name, fmt in SPC_HEADER_FIELDS:
            size = struct.calcsize(fmt)
            values = struct.unpack_from(struct_prefix + fmt, buffer, offset)
            value = values[0] if len(values) == 1 else values

            if isinstance(value, bytes):
                value = value.split(b"\x00")[0].decode("latin-1", errors="replace").strip()

            header[field_name] = value
            offset += size

        return header

    def _read_x_axis(self, f) -> np.ndarray:
        """Read or generate X coordinate array."""
        flags = self.header["flags"]
        n_points = self.header["n_points"]

        # Skip to position after main header
        f.seek(SPC_HEADER_SIZE)

        if flags & FLAG_EXPLICIT_X and not (flags & FLAG_PER_SUBFILE_XY):
            # Explicit global X array, this is stored directly after main header
            dt = "<f4" if self._struct_prefix == "<" else ">f4"
            byte_count = int(n_points) * 4
            buf = f.read(byte_count)
            if len(buf) != byte_count:
                raise ValueError(
                    f"Could not read explicit global X axis: expected {byte_count} bytes, got {len(buf)}"
                )
            x_data = np.frombuffer(buf, dtype=dt).astype(np.float64)
        else:
            # Evenly spaced X defined from first_x and last_x in header
            x_data = np.linspace(self.header["first_x"], self.header["last_x"], n_points)

        return x_data

    def _read_all_subfiles(self, f) -> tuple[np.ndarray, list[dict[str, object]]]:
        """Read all Y data and subheaders.

        Returns:
            y: 1D array for single spectrum, 2D [n_points, n_subfiles] for multifile
            subheaders: List of subheader dicts
        """
        flags = self.header["flags"]
        n_points = self.header["n_points"]
        n_subfiles = self.header["n_subfiles"]
        is_16bit_y = bool(flags & FLAG_Y_16BIT)

        # Position file pointer after main header and optional explicit global X
        f.seek(SPC_HEADER_SIZE)
        if flags & FLAG_EXPLICIT_X and not (flags & FLAG_PER_SUBFILE_XY):
            f.seek(n_points * 4, 1)  # Skip explicit global X array

        # Read all subfiles
        subheaders = []
        y_data_list = []

        for _ in range(n_subfiles):
            subheader = self._read_subheader(f)
            subheaders.append(subheader)

            # Determine Y exponent
            if flags & FLAG_IS_MULTIFILE:
                y_exponent = subheader["exponent"]
            else:
                y_exponent = self.header["exponent"]

            # Read Y data
            y_values = self._read_y_data(f, n_points, y_exponent, is_16bit_y)
            y_data_list.append(y_values)

        # Stack Y data
        if n_subfiles == 1:
            y_array = y_data_list[0]  # 1D for single spectrum
        else:
            y_array = np.column_stack(y_data_list)  # 2D [n_points, n_subfiles]

        return y_array, subheaders

    @staticmethod
    def _build_subfiles_shared_x(x: np.ndarray, y: np.ndarray, subheaders: list[dict[str, object]]) -> list[SPCSubfile]:
        """Build SPCSubfile instances for shared-X files."""
        if y.ndim == 1:
            # Single subfile
            return [SPCSubfile(x=x, y=y, z=subheaders[0].get("z_value"), subheader=subheaders[0])]

        subfiles: list[SPCSubfile] = []
        for i, subheader in enumerate(subheaders):
            # Multiple subfiles, extract column i from Y array
            subfiles.append(SPCSubfile(x=x, y=y[:, i], z=subheader.get("z_value"), subheader=subheader))
        return subfiles

    def _read_all_subfiles_xyxys(self, f) -> tuple[list[SPCSubfile], list[dict[str, object]]]:
        """Read XYXY (TXYXYS) files where each subfile has its own X array and length."""
        if not self.explicit_x:
            raise ValueError("XYXY (TXYXYS) files require explicit X values (TXVALS flag)")

        flags = self.header["flags"]
        n_subfiles = self.header["n_subfiles"]
        is_16bit_y = flags & FLAG_Y_16BIT

        f.seek(SPC_HEADER_SIZE)

        subheaders: list[dict[str, object]] = []
        subfiles: list[SPCSubfile] = []

        for _ in range(n_subfiles):
            subheader = self._read_subheader(f)
            subheaders.append(subheader)

            n_points = int(subheader["n_points"])
            dt_x = "<f4" if self._struct_prefix == "<" else ">f4"
            x = np.frombuffer(f.read(n_points * 4), dtype=dt_x).astype(np.float64)

            if flags & FLAG_IS_MULTIFILE:
                y_exponent = int(subheader["exponent"])
            else:
                y_exponent = int(self.header["exponent"])

            y = self._read_y_data(f, n_points, y_exponent, is_16bit_y)
            subfiles.append(SPCSubfile(x=x, y=y, z=subheader.get("z_value"), subheader=subheader))

        return subfiles, subheaders

    def _read_subheader(self, f) -> dict[str, object]:
        """Read and parse a 32-byte subheader."""
        buffer = f.read(SPC_SUBHEADER_SIZE)
        if len(buffer) < SPC_SUBHEADER_SIZE:
            raise ValueError(f"Could not read subheader: got {len(buffer)} bytes, expected {SPC_SUBHEADER_SIZE}")

        subheader: dict[str, object] = {}
        offset = 0

        for field_name, fmt in SPC_SUBHEADER_FIELDS:
            size = struct.calcsize(fmt)
            values = struct.unpack_from(self._struct_prefix + fmt, buffer, offset)
            value = values[0] if len(values) == 1 else values
            if isinstance(value, bytes):
                value = value.decode("latin-1", errors="replace").strip()
            subheader[field_name] = value
            offset += size

        return subheader

    def _read_y_data(self, f, n_points: int, exponent: int, is_16bit: bool) -> np.ndarray:
        """Read and decode Y values for one subfile."""
        # Read raw data
        if is_16bit:
            dt = "<i2" if self._struct_prefix == "<" else ">i2"
            byte_count = int(n_points) * 2
            buf = f.read(byte_count)
            if len(buf) != byte_count:
                raise ValueError(
                    f"File ended while reading Y data: expected {byte_count} bytes, got {len(buf)}"
                )
            y_raw = np.frombuffer(buf, dtype=dt)
        else:
            dt = "<i4" if self._struct_prefix == "<" else ">i4"
            byte_count = int(n_points) * 4
            buf = f.read(byte_count)
            if len(buf) != byte_count:
                raise ValueError(
                    f"File ended while reading Y data: expected {byte_count} bytes, got {len(buf)}"
                )
            y_raw = np.frombuffer(buf, dtype=dt)

        # Decode to float
        if exponent == -128:  # 0x80 = floating point
            if is_16bit:
                # Check for invalid combination
                raise ValueError("Cannot have 16-bit Y with float exponent")
            ft = "<f4" if self._struct_prefix == "<" else ">f4"
            y_data = y_raw.view(ft).astype(np.float64)
        else:
            # Fixed-point conversion
            bit_width = 16 if is_16bit else 32
            scale = 2.0**exponent / (2.0**bit_width)
            y_data = y_raw.astype(np.float64) * scale

        return y_data

    def _read_log_block(self, f) -> str | None:
        """Read the log block if present."""
        log_offset = self.header.get("log_offset", 0)
        if log_offset == 0:
            return None

        # Read LOGSTC header (64 bytes)
        f.seek(log_offset)
        logstc_buffer = f.read(64)
        if len(logstc_buffer) < 64:
            return None

        # Parse LOGSTC structure (5 integers)
        fmt = self._struct_prefix + "IIIII"
        logsize, _, txt_offset, binary_size, logdsks = struct.unpack(fmt, logstc_buffer[:20])
        # Remaining 44 bytes are reserved/spare

        # Read remaining log block data (we already read first 64 bytes)
        remaining_data = f.read(logsize - 64)
        log_data = logstc_buffer + remaining_data

        # Extract text starting at text offset
        if txt_offset >= len(log_data):
            return None

        text_data = log_data[txt_offset:]

        # Decode log text up to the first null byte
        # Log text is ASCII with CR+LF line endings, terminated by \0
        log_text = text_data.split(b"\x00")[0].decode("latin-1", errors="replace")
        return log_text
