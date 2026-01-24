# SPCpy - An SPC file reader in Python
A modern reader for GRAMS/Thermo-Galactic [SPC files](https://en.wikipedia.org/wiki/SPC_file_format) - a widely used file format in spectroscopy.


## Features
SPCpy focuses on a small, practical API for loading SPC files into NumPy arrays. It aims to make common spectroscopy workflows easy (load → inspect metadata → work with `x`/`y` arrays), while keeping the implementation straightforward and well-tested.


## Supported SPC formats
Currently supported:

- New-format SPC (512-byte header): `0x4B` (little-endian, tested) and `0x4C` (big-endian, supported but not currently tested due to lack of sample files)
- X modes: implicit evenly-spaced X, explicit global X (`TXVALS`), and per-subfile X/Y (`TXYXYS` / "XYXY")
- Multifile (`TMULTI`) and minimal 4D metadata: `.z` / `.w` expose per-subfile coordinates and `.w_planes` exposes the plane count
- Log text block is read when present

## Usage
```python
from spcpy import SPCFile

# Load a file
spc = SPCFile("path/to/your/file.spc")

# Single-spectrum (or single-subfile XYXY): x and y are 1D
x = spc.x
y = spc.y

# Multifile with shared X: y is (n_points, n_subfiles)
if spc.is_multifile and spc.has_shared_x:
    y0 = spc.y[:, 0]

# Iterate subfiles (always works)
for sub in spc:
    print(sub.z, sub.x.shape, sub.y.shape)

# Per-subfile XYXY (TXYXYS) multifiles: access per-spectrum arrays
if spc.per_subfile_xy and spc.is_multifile:
    x0 = spc[0].x
    y0 = spc[0].y
```

Log text (if present):

```python
spc = SPCFile("path/to/ftir.spc")
if spc.log:
    print(spc.log.splitlines()[0])
```

## Limitations
SPCpy currently rejects old-format `0x4D` SPC files. For `TXYXYS` files, SSFSTC directory-based random access is not implemented yet.

