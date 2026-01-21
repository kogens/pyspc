# PySPC - An SPC file reader in Python
A modern reader for GRAMS/Thermo-Galactic [SPC files](https://en.wikipedia.org/wiki/SPC_file_format).

## Features
- Read single and multi-file SPC files (1D, 2D, 3D spectra)
- Access spectral data, axes, and metadata
- Support for various SPC subfile types (FTIR, Raman, UV-Vis, etc.)
- Easy integration with NumPy for data analysis
- Comprehensive test suite using pytest
- Open source under the LGPL-3.0 License


## Usage
```python
from pyspc import SPCFile
import matplotlib.pyplot as plt

# Load an SPC file (shared x-axis, single or multifile)
spc = SPCFile("path/to/your/file.spc")
plt.plot(spc.x, spc.y)

# Indexing into multifile spectra
subfile = spc[0]
print(subfile)
```




