# SPCfile - An SPC file reader in Python
A modern reader for GRAMS/Thermo-Galactic [SPC files](https://en.wikipedia.org/wiki/SPC_file_format) - a widely used file format in spectroscopy.

## Installation
Requires Python >= 3.9 and NumPy. Install directly from this repository using `pip`:

```bash
pip install git+https://github.com/kogens/spcfile.git
```


## Features
SPCfile focuses on a small, practical API for loading SPC files into NumPy arrays. It aims to make common spectroscopy workflows easy (load → inspect metadata → work with `x`/`y` arrays), while keeping the implementation straightforward and well-tested.


## Supported SPC formats
Currently supported:

- New-format SPC (512-byte header): `0x4B` (little-endian, tested) and `0x4C` (big-endian, supported but not currently tested due to lack of sample files)
- X modes: implicit evenly-spaced X, explicit global X (`TXVALS`), and per-subfile X/Y (`TXYXYS` / "XYXY")
- Multifile (`TMULTI`) and minimal 4D metadata: `.z` / `.w` expose per-subfile coordinates and `.w_planes` exposes the plane count
- Log text block is read when present

## Usage examples
```python
from spcfile import SPCFile

# Load multifile with a shared X axis
spc = SPCFile("multifile.spc")
print(spc)
```

```
SPC File: multifile.spc
Date: 2021-03-15 09:21:00
Subfiles: 1
Points per subfile: 1776
Experiment type: General
Units: X='Wavenumber (cm^-1)', Y='Transmittance', Z='Arbitrary', W='Arbitrary'
```

Datapoints and coordinates are easily accessible as attributes of the `SPCFile` object:

```python
# Access X, Y, Z values (represented as numpy arrays)
x = spc.x
y = spc.y
z = spc.z

# Subfiles can be indexed like a list
subfile0 = spc[0]
y0 = subfile0.y
```


## Limitations
SPCfile currently rejects old-format `0x4D` SPC files.

