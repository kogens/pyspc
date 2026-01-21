"""GRAMS SPC file reader."""

from __future__ import annotations

import struct
from dataclasses import dataclass
import datetime
from pathlib import Path
from typing import Iterable

import numpy as np

# Main header flag bits (ftflgs)
FLAG_Y_16BIT = 0x01  # Y data is 16-bit if set, 32-bit if clear
FLAG_CHROMATOGRAM = 0x02  # Enables chromatogram / fexper interpretation
FLAG_MULTIFILE = 0x04  # Multifile; more than one subfile present
FLAG_RANDOM_Z = 0x08  # Multifile with arbitrary (unordered) Z values
FLAG_ORDERED_Z = 0x10  # Multifile with ordered but uneven Z values
FLAG_CUSTOM_AXIS_LABELS = 0x20  # Use fcatxt axis labels instead of type defaults
FLAG_PER_SUBFILE_XY = 0x40  # Per-subfile X arrays and lengths (TXYXYS mode)
FLAG_EXPLICIT_X = 0x80  # X values stored explicitly as float array(s)


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

SPC_HEADER_SIZE = 512
SPC_SUBHEADER_SIZE = 32

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
            self.header = self._read_header(f)

            if self.header["version"] != 0x4B:
                raise ValueError(f"Unsupported SPC version: {self.header['version']:02X}")

            if self.header.get("w_planes", 0) != 0:
                raise NotImplementedError("4D W-plane data is not supported yet")

            # Read X axis
            self._x = self._read_x_axis(f)

            # Read all Y data and subheaders
            self._y, self.subheaders = self._read_all_subfiles(f)

            # Read log block
            self.log = self._read_log_block(f)

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
        hour   = (dt_raw >> 6) & 0x1F
        day    = (dt_raw >> 11) & 0x1F
        month  = (dt_raw >> 16) & 0x0F
        year   = (dt_raw >> 20) & 0xFFF

        # Basic validation
        if not ((1 <= month <= 12) and (1 <= day <= 31) and (1900 < year < 2100)):
            return None

        return datetime.datetime(year, month, day, hour, minute)

    @property
    def has_shared_x(self) -> bool:
        """True if all subfiles share a common X axis.

        False for TXYXYS files where each subfile has its own X array.
        """
        flags = self.header["flags"]
        return not (flags & FLAG_PER_SUBFILE_XY)

    @property
    def x(self) -> np.ndarray:
        """X coordinates (only for files with shared X axis).

        Returns:
            1D array of X coordinates shared across all subfiles.

        Raises:
            ValueError: For TXYXYS files with per-subfile X arrays.
                        Use spc[i].x to access individual X arrays.
        """
        if not self.has_shared_x:
            raise ValueError(
                "This SPC file has per-subfile X arrays (TXYXYS mode). Use spc[i].x to access X for each spectrum."
            )
        return self._x

    @property
    def y(self) -> np.ndarray:
        """Y values (only for files with shared X axis).

        Returns:
            1D array for single spectrum.
            2D array [n_points, n_subfiles] for multifile.

        Raises:
            ValueError: For TXYXYS files with varying lengths.
                        Use spc[i].y to access individual Y arrays.
        """
        if not self.has_shared_x:
            raise ValueError(
                "This SPC file has per-subfile XY arrays (TXYXYS mode). Use spc[i].y to access Y for each spectrum."
            )
        return self._y

    def __len__(self) -> int:
        """Number of subfiles (spectra) in the file."""
        return len(self.subheaders)

    def __getitem__(self, index: int) -> SPCSubfile:
        """Get k-th spectrum as an SPCSubfile."""
        y_data = self._y[:, index] if self._y.ndim == 2 else self._y
        z_value = self.subheaders[index].get("z_value")
        return SPCSubfile(x=self._x, y=y_data, z=z_value, subheader=self.subheaders[index])

    def __iter__(self) -> Iterable[SPCSubfile]:
        """Iterate over all subfiles."""
        for i in range(len(self)):
            yield self[i]

    def __repr__(self) -> str:
        flags = self.header["flags"]
        parts = ["multifile" if flags & FLAG_MULTIFILE else "single"]
        parts.append("shared-x" if self.has_shared_x else "per-subfile-x")
        if flags & FLAG_EXPLICIT_X:
            parts.append("explicit-x")
        flags = "|".join(parts)

        return f"<SPCFile path={self._path!r} n_subfiles={len(self)} n_points={self.header['n_points']} flags={flags}>"

    def __str__(self) -> str:
        lines = [
            f"SPC File: {self._path}",
            f"Date: {self.date}",
            f"Subfiles: {len(self)}",
            f"Points per subfile: {self.header['n_points']}",
            f"Experiment type: {self.experiment}",
            f"X unit: {self.x_unit}",
            f"Y unit: {self.y_unit}",
            f"Z unit: {self.z_unit}",
            f"W unit: {self.w_unit}",
        ]
        return "\n".join(lines)

    def _read_header(self, f) -> dict[str, object]:
        """Read and parse the 512-byte main header."""
        f.seek(0)
        buffer = f.read(SPC_HEADER_SIZE)
        if len(buffer) < SPC_HEADER_SIZE:
            raise ValueError(f"File too small: got {len(buffer)} bytes, expected {SPC_HEADER_SIZE}")

        header: dict[str, object] = {}
        offset = 0

        for field_name, fmt in SPC_HEADER_FIELDS:
            size = struct.calcsize(fmt)
            values = struct.unpack_from(fmt, buffer, offset)
            value = values[0] if len(values) == 1 else values

            # Decode byte strings
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
            x_data = np.frombuffer(f.read(n_points * 4), dtype="<f4").astype(np.float64)
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
            if flags & FLAG_MULTIFILE:
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

    def _read_subheader(self, f) -> dict[str, object]:
        """Read and parse a 32-byte subheader."""
        buffer = f.read(SPC_SUBHEADER_SIZE)
        if len(buffer) < SPC_SUBHEADER_SIZE:
            raise ValueError(f"Could not read subheader: got {len(buffer)} bytes, expected {SPC_SUBHEADER_SIZE}")

        subheader: dict[str, object] = {}
        offset = 0

        for field_name, fmt in SPC_SUBHEADER_FIELDS:
            size = struct.calcsize(fmt)
            values = struct.unpack_from(fmt, buffer, offset)
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
            y_raw = np.frombuffer(f.read(n_points * 2), dtype="<i2")
        else:
            y_raw = np.frombuffer(f.read(n_points * 4), dtype="<i4")

        # Decode to float
        if exponent == -128:  # 0x80 = floating point
            if is_16bit:
                # Check for invalid combination
                raise ValueError("Cannot have 16-bit Y with float exponent")
            y_data = y_raw.view("<f4").astype(np.float64)
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
        logsize, _, txt_offset, binary_size, logdsks = struct.unpack("<IIIII", logstc_buffer[:20])
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
