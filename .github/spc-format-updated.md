# GRAMS SPC File Format

This document summarizes the GRAMS SPC binary format as exposed by Galactic's SPC SDK (SPC.H) and Universal Data File specification. It is aimed at developers implementing SPC readers/writers, with emphasis on single-spectrum and simple multifile XY data.

> The canonical reference remains the original SPC SDK headers and UDF specification. This document is a practical overview, not a verbatim copy.

---

## 1. High-level overview

An SPC file is a binary container for one or more spectra ("subfiles"):

1. **Main header**: fixed 512-byte structure (SPCHDR) in the "new" format.
2. **Optional global X array**: `fnpts` 32-bit floats if TXVALS is set and TXYXYS is clear.
3. **One or more subfiles**:
   - SUBHDR (32-byte per subfile header).
   - Optional per-subfile X array (when TXYXYS is set).
   - Y data array (16-bit or 32-bit fixed-point or 32-bit IEEE float).
4. **Optional XY directory** for TXYXYS multifiles (array of SSFSTC).
5. **Optional log block** containing a LOGSTC header, optional binary data, and log text.

There is also an older "0x4D" format with a shorter header (224 or 256 bytes) and some different field types. Modern software normally writes the new 512-byte header form.

### File naming conventions

By convention, SPC files use specific extensions to indicate data type:

| Extension | Data Type |
|-----------|-----------|
| `.SPC` | Spectrum (general case) |
| `.CGM` | Chromatogram (general case) |
| `.GC` | Gas Chromatogram |
| `.LC` | Liquid Chromatogram |
| `.HC` | HPLC Chromatogram |
| `.FIR` | Infrared (FT-IR) Spectrum |
| `.IR` | Near Infrared (NIR) Spectrum |
| `.VIS` | Visible Spectrum |
| `.UV` | Ultraviolet (Visible) Spectrum |
| `.XRY` | X-ray Spectrum |
| `.MS` | Mass Spectrum |
| `.NMR` | Nuclear Magnetic Resonance Spectrum |
| `.RMN` | Raman Spectrum |
| `.EFL` | Fluorescence Spectrum |
| `.AAS` | Atomic Spectrum |
| `.DAS` | Diode Array Spectrum |

For hyphenated techniques (e.g., GC-IR, GC-MS), GRAMS software automatically "links" files if spectral data is in a `.SPC` file and chromatogram data is in a `.CGM` file with the same base name in the same directory.

---

## 2. Endianness and versions

The main header includes a version/endianness byte `fversn`:

- `0x4B`: new format, little-endian (LSB first). This is the common case.
- `0x4C`: new format, big-endian (MSB first).
- `0x4D`: old format (header 224/256 bytes, some fields differ, 32-bit Y words are word-swapped).

For a Python parser using `struct`, you typically:

- Read the 512-byte header using little-endian (`<`) and check that `fversn == 0x4B` for now.
- Optionally add support for `0x4D` later (see section 13).

The new format supports multifiles, the Audit Log, and has a 512-byte main header. The old format does not support these features and uses a 256-byte main header.

---

## 3. Basic C types and Python `struct` formats

SPC structures use fixed-width C types with no padding between fields. Common types and their sizes:

| C type      | Size (bytes) | Python `struct` (little-endian) |
| ----------- | ------------ | -------------------------------- |
| `BYTE`      | 1            | `B`                              |
| `char`      | 1            | `c` or `s` (for arrays)          |
| `WORD`      | 2 (uint16)   | `H`                              |
| `DWORD`     | 4 (uint32)   | `I`                              |
| `float`     | 4 (IEEE-754) | `f`                              |
| `double`    | 8 (IEEE-754) | `d`                              |
| `char[n]`   | n            | `ns` (e.g. `9s`, `130s`)         |

Use a leading `<` or `>` in the format string depending on endianness.

Example for reading a `DWORD` followed by a `double` in Python:

```python
npts, first_x = struct.unpack('<Id', f.read(12))
```

---

## 4. Main header (SPCHDR)

The main header is a packed 512-byte structure with this C layout (simplified names/ordering only):

```c
typedef struct {
    BYTE   ftflgs;       // flags (TSPREC, TMULTI, TXVALS, ...)
    BYTE   fversn;       // 0x4B=new LSB, 0x4C=new MSB, 0x4D=old
    BYTE   fexper;       // instrument / technique code
    char   fexp;         // Y exponent; 0x80 (-128) => float Y
    DWORD  fnpts;        // # of points, or XY directory offset for TXYXYS
    double ffirst;       // X of first point
    double flast;        // X of last point
    DWORD  fnsub;        // # of subfiles (1 if not TMULTI)
    BYTE   fxtype;       // X axis unit type
    BYTE   fytype;       // Y axis unit type
    BYTE   fztype;       // Z axis unit type
    BYTE   fpost;        // posting disposition
    DWORD  fdate;        // packed date/time
    char   fres[9];      // resolution string
    char   fsource[9];   // instrument string
    WORD   fpeakpt;      // peak point for interferograms (0 if unknown)
    float  fspare[8];    // reserved / internal (used by Array Basic)
    char   fcmnt[130];   // comment text
    char   fcatxt[30];   // axis labels if TALABS
    DWORD  flogoff;      // offset of log block, or 0
    DWORD  fmods;        // modification flags
    BYTE   fprocs;       // processing code
    BYTE   flevel;       // calibration level + 1
    WORD   fsampin;      // sample injection number
    float  ffactor;      // multiplier / concentration factor
    char   fmethod[48];  // method / program / data filenames
    float  fzinc;        // Z increment (0 => use first subfile delta)
    DWORD  fwplanes;     // # of W-planes (4D data), or 0
    float  fwinc;        // W increment (if fwplanes != 0)
    BYTE   fwtype;       // W axis unit type
    char   freserv[187]; // reserved (must be zero)
} SPCHDR;
```

Field summary:

| Name      | Type        | `struct` code | Description                         |
| --------- | ----------- | ------------- | ----------------------------------- |
| `ftflgs`  | `uint8`     | `B`           | Flags (precision, multifile, X mode) |
| `fversn`  | `uint8`     | `B`           | Format / endianness (0x4B, 0x4C, 0x4D) |
| `fexper`  | `uint8`     | `B`           | Instrument / technique code          |
| `fexp`    | `int8`      | `b`           | Y exponent (0x80 => float Y values)  |
| `fnpts`   | `uint32`    | `I`           | Point count (or XY dir offset TXYXYS) |
| `ffirst`  | `double`    | `d`           | X coordinate of first point          |
| `flast`   | `double`    | `d`           | X coordinate of last point           |
| `fnsub`   | `uint32`    | `I`           | Number of subfiles                   |
| `fxtype`  | `uint8`     | `B`           | X axis unit type                     |
| `fytype`  | `uint8`     | `B`           | Y axis unit type                     |
| `fztype`  | `uint8`     | `B`           | Z axis unit type                     |
| `fpost`   | `uint8`     | `B`           | Posting disposition (see 4.2.1)      |
| `fdate`   | `uint32`    | `I`           | Packed date/time                     |
| `fres`    | `char[9]`   | `9s`          | Resolution text                      |
| `fsource` | `char[9]`   | `9s`          | Source instrument text               |
| `fpeakpt` | `uint16`    | `H`           | Peak point index (interferograms)    |
| `fspare`  | `float[8]`  | `8f`          | Reserved / internal                  |
| `fcmnt`   | `char[130]` | `130s`        | Comment text                         |
| `fcatxt`  | `char[30]`  | `30s`         | Axis label strings (if TALABS)       |
| `flogoff` | `uint32`    | `I`           | Offset of log block                  |
| `fmods`   | `uint32`    | `I`           | Modification flags                   |
| `fprocs`  | `uint8`     | `B`           | Processing code                      |
| `flevel`  | `uint8`     | `B`           | Calibration level + 1                |
| `fsampin` | `uint16`    | `H`           | Sample injection number              |
| `ffactor` | `float`     | `f`           | Data multiplier / concentration      |
| `fmethod` | `char[48]`  | `48s`         | Method/program filenames             |
| `fzinc`   | `float`     | `f`           | Z increment                          |
| `fwplanes`| `uint32`    | `I`           | Number of W-planes (4D)              |
| `fwinc`   | `float`     | `f`           | W increment                          |
| `fwtype`  | `uint8`     | `B`           | W axis unit type                     |
| `freserv` | `char[187]` | `187s`        | Reserved (zeros)                     |

### 4.1 Python `struct` format for SPCHDR (little-endian)

A minimal format string for the fields above in Python (omitting any padding) is:

```python
SPCHDR_FMT = '<BBBBI d d I BBBB I 9s 9s H 8f 130s 30s I I B B H f 48s f I f B 187s'
```

You can unpack with:

```python
hdr = struct.unpack(SPCHDR_FMT, f.read(512))
```

In practice you will want a small helper that maps the tuple into a named object or `dataclass`.

### 4.2 Key header fields

For a reader, the most important SPCHDR fields are:

- `ftflgs` (BYTE): feature flags bitfield:

   | Flag     | Hex   | Meaning                                                |
   | -------- | ----- | ------------------------------------------------------ |
   | `TSPREC` | 0x01  | Y data is 16-bit if set, 32-bit if clear               |
   | `TCGRAM` | 0x02  | Enables chromatogram / `fexper` interpretation        |
   | `TMULTI` | 0x04  | Multifile; more than one subfile present              |
   | `TRANDM` | 0x08  | Multifile with arbitrary (unordered) Z values         |
   | `TORDRD` | 0x10  | Multifile with ordered but uneven Z values            |
   | `TALABS` | 0x20  | Use `fcatxt` axis labels instead of type defaults     |
   | `TXYXYS` | 0x40  | Per-subfile X arrays and lengths (TXYXYS mode)        |
   | `TXVALS` | 0x80  | X values stored explicitly as float array(s)          |

- `fversn` (BYTE): version and endianness; usually `0x4B`.
- `fexper` (BYTE): instrument / technique (FT-IR, MS, NMR, etc.).
- `fexp` (char, signed): Y scaling exponent (see section 6).
- `fnpts` (DWORD):
   - Normal case (no TXYXYS): number of points per subfile.
   - TXYXYS multifile: byte offset of the XY directory (array of SSFSTC).
- `ffirst`, `flast` (double): X coordinate range of the first subfile.
- `fnsub` (DWORD): number of subfiles; 1 for single-spectrum files.
- `fxtype`, `fytype`, `fztype` (BYTE): axis unit types (see 4.3).
- `flogoff` (DWORD): offset of the optional log block.
- `fwplanes`, `fwinc`, `fwtype`: 4D (W-axis) information (see 11).

Other fields (`fspare`, `fmethod`, etc.) can often be safely ignored for basic reading.

**Important notes:**

- For XY files (TXVALS set), `ffirst` and `flast` do not necessarily correspond to the actual X-axis limits. You should get the actual limits by reading the X array values.
- For TXYXYS multifiles, you must scan through all subfiles and get the outermost X limits, handling both low-to-high and high-to-low data.
- `fpeakpt` is the data point index (0-based) of the ZPD point for interferograms. Set to 0 if not known or not an interferogram.

#### 4.2.1 Posting disposition (`fpost`)

The `fpost` field specifies the desired post-collection processing and file storage behavior for the data in GRAMS software. Zero values indicate unspecified/default settings. The field uses constants defined in GRAMSDDE.H:

| Value | Constant  | Meaning |
| ----- | --------- | ------- |
| 0     | `PSTDEFT` | Use default setting (unspecified) |
| 1     | `PSTSAVE` | Save file to disk (but remove from memory) |
| 2     | `PSTAPPD` | Append to end of database |
| 3     | `PSTMERG` | Merge into current database row |
| 4     | `PSTBACK` | Save as new background for experiment |
| 5     | `PSTNONE` | Do not save after processing |
| 6     | `PSTKEEP` | Do not save and keep in memory as #S |
| 7     | `PSTBOTH` | Both disk save & keep in memory (ABC Driver Only) |
| 8     | `PSTASK`  | Query user: save, keep, or both (ABC Driver Only) |

> **Note for readers**: This field was primarily used by GRAMS software for workflow control and does not affect how spectral data is parsed or interpreted. Modern readers can treat it as metadata only.

#### 4.2.2 Packed date format (`fdate`)

The date/time is encoded as unsigned integers into a 32-bit value (most significant bit on left):

- **Year**: 12 bits (0-4095)
- **Month**: 4 bits (1-12)
- **Day**: 5 bits (1-31)
- **Hour**: 5 bits (0-23)
- **Minute**: 6 bits (0-59)

Example Python decoding:

```python
minute = fdate & 0x3F
hour = (fdate >> 6) & 0x1F
day = (fdate >> 11) & 0x1F
month = (fdate >> 16) & 0x0F
year = (fdate >> 20) & 0xFFF
```

### 4.3 Axis unit types (X, Y, Z, W)

The headers define enumerations for axis types. Common values:

**X axis** (`fxtype`, `fztype`, `fwtype`):

| Code | Meaning                 |
| ---- | ----------------------- |
| 0    | Arbitrary               |
| 1    | Wavenumber (cm⁻¹)       |
| 2    | Micrometers             |
| 3    | Nanometers              |
| 4    | Seconds                 |
| 5    | Minutes                 |
| 30   | Hours                   |
| 6    | Hertz (Hz)              |
| 7    | Kilohertz (kHz)         |
| 8    | Megahertz (MHz)         |
| 26   | Gigahertz (GHz)         |
| 9    | Mass (m/z)              |
| 13   | Raman shift (cm⁻¹)      |
| 22   | Data points             |
| 18   | Degrees                 |
| 19   | Temperature (F)         |
| 20   | Temperature (C)         |
| 21   | Temperature (K)         |
| 27   | Centimeters (cm)        |
| 28   | Meters (m)              |
| 29   | Millimeters (mm)        |
| 23   | Milliseconds (mSec)     |
| 24   | Microseconds (uSec)     |
| 25   | Nanoseconds (nSec)      |

**Y axis** (`fytype`):

| Code | Meaning                 |
| ---- | ----------------------- |
| 0    | Arbitrary intensity     |
| 1    | Interferogram           |
| 2    | Absorbance              |
| 3    | Kubelka-Munk            |
| 10   | Log(1/R)                |
| 11   | Percent                 |
| 4    | Counts                  |
| 5    | Volts                   |
| 9    | Millivolts (mV)         |
| 12   | Intensity               |
| 13   | Relative intensity      |
| 128  | Transmission (valleys)  |
| 129  | Reflectance             |
| 130  | Arbitrary or Single Beam with Valley Peaks |
| 131  | Emission                |

**Note**: All Y types ≥128 are assumed to have inverted (valley) peaks.

For full lists and rare values, see the original header definitions.

#### Custom axis labels

If the `TALABS` flag is set in `ftflgs`, axis labels are taken from the `fcatxt` field instead of the enumerated types. The `fcatxt` contains three null-terminated strings for X, Y, and Z labels in that order. Each label can be up to 20 characters, and all three must fit within 30 bytes total.

If a label in `fcatxt` is just a null byte (empty string), the corresponding enumerated type label is used instead.

---

## 5. Subfile header (SUBHDR)

Every file, even single-spectrum files, has at least one subfile header immediately before the Y data for that subfile. In C:

```c
typedef struct {
    BYTE  subflgs;   // changed / no-peak-table / modified flags
    char  subexp;    // Y exponent for this subfile (0x80 => float)
    WORD  subindx;   // subfile index (0 = first)
    float subtime;   // Z coordinate (time, etc.) for this subfile
    float subnext;   // Z coordinate for next subfile
    float subnois;   // noise estimate (high byte non-zero if valid)
    DWORD subnpts;   // # points in TXYXYS mode; ignored otherwise
    DWORD subscan;   // # co-added scans or 0
    float subwlevel; // W-axis value for this subfile (4D)
    char  subresv[4];// reserved (zero)
} SUBHDR;
```

Field summary:

| Name        | Type        | `struct` code | Description                          |
| ----------- | ----------- | ------------- | ------------------------------------ |
| `subflgs`   | `uint8`     | `B`           | Subfile flags (changed, etc.)        |
| `subexp`    | `int8`      | `b`           | Y exponent for this subfile          |
| `subindx`   | `uint16`    | `H`           | Subfile index (0 = first)            |
| `subtime`   | `float`     | `f`           | Z coordinate (time, etc.)            |
| `subnext`   | `float`     | `f`           | Z coordinate of next subfile         |
| `subnois`   | `float`     | `f`           | Noise estimate                       |
| `subnpts`   | `uint32`    | `I`           | Point count in TXYXYS mode           |
| `subscan`   | `uint32`    | `I`           | Number of co-added scans             |
| `subwlevel` | `float`     | `f`           | W-axis value for this subfile        |
| `subresv`   | `char[4]`   | `4s`          | Reserved (zeros)                     |

Python `struct` format (little-endian):

```python
SUBHDR_FMT = '<B b H f f f I I f 4s'
```

Key fields for readers:

- `subexp`: Y exponent for this subfile. If `TMULTI` is clear, `fexp` usually applies instead.
- `subindx`: Must be correct for all subfiles, starting from 0.
- `subtime`: Z-axis coordinate (e.g. time, spectrum index).
- `subnext`: Ending Z coordinate for this subfile. Often equals `subtime`, but can differ if data collection spans a time range.
- `subnois`: Noise level used by GRAMS peak picking. The high byte must be non-zero for this value to be considered valid.
- `subnpts`: number of XY points **only** in TXYXYS multifiles.
- `subscan`: number of scans used (can be exposed as metadata).
- `subwlevel`: W-axis value for 4D data (when `fwplanes != 0`).

`subflgs` has a few bits defined:

| Bit | Hex  | Meaning |
|-----|------|---------|
| 1   | 0x01 | Subfile changed |
| 8   | 0x08 | Do not use peak table file |
| 128 | 0x80 | Subfile modified by arithmetic |

### Important notes on Z values

For evenly-spaced multifiles (`TMULTI` set, but `TORDRD` and `TRANDM` clear):

- `subtime` and `subnext` are optional (and ignored) for all but the first subfile.
- The Z spacing is determined by either:
  - The difference `subnext - subtime` from the first subfile, OR
  - The `fzinc` value in the main header if non-zero.
- If `fzinc` is zero, you must use the first subfile's `subnext - subtime` to calculate spacing.

For ordered or random multifiles (`TORDRD` or `TRANDM` set):

- Every subfile must have a meaningful `subtime` value.
- `subnext` is typically set equal to `subtime`.

---

## 6. Y data encoding

Y values are stored as fixed-point fractions with a shared exponent, or as native floats.

### 6.1 Fixed-point representation

For most files, Y data is fixed-point signed values scaled by a power of 2. Let:

- `exp` be the signed integer exponent (`fexp` from SPCHDR, or `subexp` for multifiles when TMULTI is set).
- `N` be the bit width of each stored Y sample (16 or 32 bits, signed two's complement).
- `ival` be the stored integer for one Y sample.

Then the corresponding floating value is:

$$y = 2^{e} \cdot \frac{i}{2^n}$$

- If `TSPREC` is **set**, Y samples are **16-bit** signed values (`int16`).
- If `TSPREC` is **clear**, Y samples are **32-bit** signed values (`int32`).

**Note**: Do NOT set `TSPREC` if `fexp == 0x80` (floating point Y data).

In Python, for 32-bit Y data with a common exponent:

```python
scale = 2.0 ** exp / (2.0 ** 32)
y = ival * scale
```

For 16-bit Y data:

```python
scale = 2.0 ** exp / (2.0 ** 16)
y = ival * scale
```

### 6.2 Floating-point Y data

If the exponent byte is `0x80` (i.e. -128 signed), then Y values on disk are stored directly as 32-bit IEEE floats, and the integer/fixed-point interpretation does not apply. This can be signaled by:

- `fexp == 0x80` for single-spectrum / uniform multifile.
- `subexp == 0x80` in SUBHDR for per-subfile exponents.

**Important**: As of December 1996, all new Galactic-written data file conversion routines output SPC files using 32-bit IEEE floating point format for Y data unless more precision is required. GRAMS and other Galactic software can read any format but newer converters always save as floating point.

Older GRAMS/LabCalc format converters always saved in fixed integer format.

### 6.3 Data layout per subfile

For a given subfile, the sequence is:

1. SUBHDR.
2. Optional X array (only when TXYXYS is set; see 10.5).
3. Y data array of length `n_pts`:
   - `n_pts = fnpts` in most cases.
   - `n_pts = subnpts` in TXYXYS multifiles.
   - Each Y sample is 2 or 4 bytes, signed, or 4-byte float when exponent is 0x80.

The subfile's Y exponent is:

- `subexp` when `TMULTI` is set and using the new format.
- `fexp` (main header) for simple single-spectrum files.

---

## 7. X values and coordinate systems

### 7.1 Evenly spaced X (no TXVALS)

When the `TXVALS` flag (bit 0x80 of `ftflgs`) is **clear**, X values are implicit and evenly spaced between the first and last X values (`ffirst` and `flast`):

- `fnpts` is the number of points per subfile.
- For point index $i \in [0, N-1]$ with $N = \text{fnpts}$ and endpoints $x_0 = \text{ffirst}$, $x_{N-1} = \text{flast}$:

$$x_i = x_0 + i \cdot \frac{x_{N-1} - x_0}{N - 1}$$

No X array is stored in the file. This is the standard format for Y-only files and is the most common format.

### 7.2 Global X array (TXVALS, no TXYXYS)

When `TXVALS` is set and `TXYXYS` is **clear**:

- Immediately after the 512-byte SPCHDR is a global X array of length `fnpts`.
- Each X value is a 32-bit IEEE float (4 bytes).
- For multifiles with this layout, all subfiles share this global X array.
- The `ffirst` and `flast` values in the main header should be set to match the first and last values in the X array.

The file layout becomes:

1. SPCHDR.
2. Global X array: `fnpts` floats.
3. For subfile 0..`fnsub-1`:
   - SUBHDR.
   - Y array of `fnpts` values.

This format is called "XYY" for multifiles (one X array, Y arrays for each subfile).

### 7.3 Per-subfile X arrays (TXYXYS)

When both `TXVALS` and `TXYXYS` are set:

- Each subfile has its own X array and can have a different number of points.
- In this mode:
  - The `subnpts` field in SUBHDR gives the number of X/Y points for that subfile.
  - The SPCHDR's `fnpts` field usually holds **the byte offset of the XY directory**, not a point count.
  - If `fnpts` is zero, no directory exists and you must scan sequentially.

Layout:

1. SPCHDR.
2. For each subfile (0..`fnsub-1`):
   - SUBHDR.
   - X array: `subnpts` floats (4-byte each).
   - Y array: `subnpts` samples (16 or 32-bit or float).
3. Optional directory (if `fnpts != 0`):
   - At byte offset `fnpts` from file start.
   - Contains `fnsub` entries of type SSFSTC.

The directory lets GRAMS locate subfiles quickly and allows subfiles to be non-contiguous in the file.

This format is called "XYXY" for multifiles. It is typically used for GC-MS data where each scan may have a different number of points.

**Important**: When `TXYXYS` is set, GRAMS software will draw data points as "sticks" (originating from Y=0) rather than connected lines, which is the traditional representation for mass spectra.

---

## 8. XY directory entries (SSFSTC)

In TXYXYS multifiles with a directory, each entry is:

```c
typedef struct {
    DWORD ssfposn;  // byte offset of SUBHDR (start of subfile)
    DWORD ssfsize;  // byte size of subfile incl. SUBHDR, X, Y
    float ssftime;  // Z time (subtime) of this subfile
} SSFSTC;
```

Field summary:

| Name      | Type     | `struct` code | Description                              |
| --------- | -------- | ------------- | ---------------------------------------- |
| `ssfposn` | `uint32` | `I`           | Byte offset of SUBHDR (start of subfile) |
| `ssfsize` | `uint32` | `I`           | Byte size of subfile (SUBHDR + X + Y)    |
| `ssftime` | `float`  | `f`           | Z time (subtime) of this subfile         |

Python format (little-endian):

```python
SSFSTC_FMT = '<I I f'
```

To read the directory:

1. Seek to offset `fnpts` from the start of the file.
2. Read `fnsub` entries using SSFSTC.
3. For each entry, you can directly seek to `ssfposn`, read a SUBHDR, and then the X/Y data.

When `fnpts == 0`, there is no directory; software can scan sequentially through the subfiles to build its own index. This allows GRAMS to add points to subfiles by moving them to the end of the file and updating the directory.

---

## 9. Log block structure (flogoff)

SPC files can contain an audit/log text block used by GRAMS for history and comments. The main header field `flogoff` is:

- `0` if no log block is present.
- Otherwise, a byte offset from the beginning of the file to a LOGSTC header.

The LOGSTC structure:

```c
typedef struct {
    DWORD logsizd;     // total size of disk log block (bytes)
    DWORD logsizm;     // size of in-memory block (>= logsizd)
    DWORD logtxto;     // offset from start of LOGSTC to text
    DWORD logbins;     // size of binary area after LOGSTC
    DWORD logdsks;     // size of disk-only area after binary
    char  logspar[44]; // reserved (zero)
} LOGSTC;
```

Field summary:

| Name       | Type        | `struct` code | Description                              |
| ---------- | ----------- | ------------- | ---------------------------------------- |
| `logsizd`  | `uint32`    | `I`           | Total size of disk log block (bytes)     |
| `logsizm`  | `uint32`    | `I`           | Size of in-memory block (>= logsizd)     |
| `logtxto`  | `uint32`    | `I`           | Offset from LOGSTC start to log text     |
| `logbins`  | `uint32`    | `I`           | Size of binary area after LOGSTC         |
| `logdsks`  | `uint32`    | `I`           | Size of disk-only area after binary      |
| `logspar`  | `char[44]`  | `44s`         | Reserved (zeros)                         |

Python format (little-endian):

```python
LOGSTC_FMT = '<I I I I I 44s'
```

### 9.1 Log block layout

The log block consists of three parts:

1. **LOGSTC header** (64 bytes)
2. **Log Binary block** (optional, `logbins` bytes)
   - Used to store binary data required for processing but not displayed
   - For NMR files, this contains the imaginary data from quadrature detection
   - This data is loaded into memory when the file is opened
3. **Log Disk block** (optional, `logdsks` bytes)
   - Used for private binary data that is NOT loaded into memory
   - This data is not carried along when the file is modified and saved back to disk
4. **Log Text block** (ASCII text)
   - At offset `logtxto` from the start of LOGSTC
   - Contains acquisition parameters and audit trail

The `logsizd` describes the total size (in bytes) of the complete log block. The `logsizm` should be `logsizd` rounded up to the nearest even multiple of 4096.

### 9.2 Log Text format

The log text:

- Is ASCII text.
- Consists of lines terminated by CR+LF ("\r\n").
- Must be terminated by a final `\0` byte after the last CR+LF.
- Uses a "key=value" format for parameters.

Each parameter is stored as:

```
KeyName = ########
```

or

```
KeyName = Text
```

where "########" represents a numerical value as ASCII text, and "Text" represents a text string. The key names are NOT case-sensitive but are conventionally capitalized.

**Important**: The first line of NMR log text should always be the `MODEL` parameter to identify the instrument/software type.

Common log text parameters:

| Key | Type | Description |
|-----|------|-------------|
| `MODEL` | Text | Instrument/software name (no spaces) |
| `BEGX` | Value | X position of first trace point |
| `ENDX` | Value | X position of last trace point |
| `NPTS` | Value | Number of trace points (0=XY pairs) |
| `BEGZ` | Value | Z position of first subfile |
| `ENDZ` | Value | Z position of last subfile |
| `NSUBS` | Value | Number of subfiles |
| `XTYPE` | Value | X axis type number (as defined in SPC.H) |
| `YTYPE` | Value | Y axis type number |
| `ZTYPE` | Value | Z axis type number |
| `RES` | Text | Collection X axis resolution (-1 if unknown) |
| `NAME` | Text | Data file name (optionally with path) |
| `MEMO` | Text | Primary comment/description |
| `USER` | Text | User's or Analyst's name |
| `SCANS` | Value | Number of co-added scans |

For a simple reader, you can:

1. Seek to `flogoff` and read LOGSTC.
2. Read `logsizd` bytes total into a buffer.
3. Take the slice starting at `logtxto`, decode as ASCII up to the first `\0` byte.

Expose this text as a single string or list of lines on your `SPCFile` object.

### 9.3 NMR-specific binary data

For NMR files, the Log Binary block is used to store the imaginary component of complex data:

- Only the "real" portion is stored in the main Y data block.
- The imaginary data is stored in the binary log block as 32-bit IEEE floats.
- Both arrays must have exactly the same number of data points.

The `logbins` field specifies the size of this binary block. For NMR files:

- The main Y data and imaginary data must be the same size.
- The imaginary data can ONLY be stored as 32-bit IEEE floats (unlike the main Y data which can be fixed-point).

---

## 10. File types: single and multifile

The combination of `TMULTI`, `TXVALS`, `TXYXYS`, `TRANDM`, and `TORDRD` describes the structure.

### Quick reference: byte-level layouts

For quick visualization, here are the common file layouts:

**Single-spectrum, evenly spaced X:**
```
[SPCHDR 512B] → [SUBHDR 32B] → [Y data: n_points × (2 or 4 bytes)]
```

**Single-spectrum, explicit X:**
```
[SPCHDR 512B] → [Global X: n_points × 4B] → [SUBHDR 32B] → [Y data]
```

**Multifile with shared X:**
```
[SPCHDR 512B] → [Optional Global X] → 
  [SUBHDR₀ 32B] → [Y₀ data] →
  [SUBHDR₁ 32B] → [Y₁ data] → ...
```

**Multifile with per-subfile X/Y (TXYXYS):**
```
[SPCHDR 512B] → 
  [SUBHDR₀ 32B] → [X₀: subnpts₀ × 4B] → [Y₀: subnpts₀ × (2 or 4 bytes)] →
  [SUBHDR₁ 32B] → [X₁: subnpts₁ × 4B] → [Y₁: subnpts₁ × (2 or 4 bytes)] → ...
  [Optional Directory]
```

Y data size per point: 2 bytes if TSPREC flag set, else 4 bytes (or 4 bytes for float if exp=0x80).

---

### 10.1 Single-spectrum, evenly spaced X (simplest case)

Typical for many spectroscopy files. This is the **standard format** for SPC files.

- Flags: `TMULTI` clear, `TXVALS` clear.
- Interpretation:
  - One subfile (but a SUBHDR is still present).
  - X axis is evenly spaced between `ffirst` and `flast` with `fnpts` points.

Layout:

1. SPCHDR (512 bytes).
2. SUBHDR (32 bytes).
3. Y data array: `fnpts` samples.

To read:

- Compute X implicitly from `ffirst`, `flast`, `fnpts`.
- Read Y as 16 or 32-bit values, apply exponent.

This is also called a "Y-only file" since X values are not explicitly stored.

### 10.2 Single-spectrum, variable X (XY data)

Use this when X spacing is not uniform (e.g., Raman spectra from CCD detectors, TOF MS).

- Flags: `TXVALS` set, `TMULTI` clear, `TXYXYS` clear.
- Interpretation:
  - One subfile.
  - A single global X array of `fnpts` float32 values.

Layout:

1. SPCHDR.
2. X array: `fnpts` floats.
3. SUBHDR.
4. Y array: `fnpts` samples.

**Important**: Although X values are fully specified by the X array, `ffirst` and `flast` must also be set to the first and last values in the X array. The X values must be in either increasing or decreasing order (not random).

### 10.3 Multifile: evenly spaced Z (Y-only)

These are data cubes in (X, Z) where each subfile is a spectrum at a different Z (e.g. time, scan number). This is the most common multifile format.

- Flags: `TMULTI` set, `TXVALS` clear.
- `TRANDM` and `TORDRD` both clear.
- Interpretation:
  - Subfiles are evenly spaced along Z.
  - All subfiles have the same evenly spaced X values.
  - Only Y values are stored for each subfile.
  - Z increment determined by either:
    - `fzinc` in main header (if non-zero), OR
    - `subnext - subtime` from first subfile.

Layout:

1. SPCHDR.
2. For each subfile `k` in 0..`fnsub-1`:
   - SUBHDR.
   - Y array: `fnpts` samples.

The initial Z value is taken from `subtime` of the first subfile. Subsequent Z values are:

$Z_n = Z_0 + (n \times \text{fzinc})$

where $n$ is the subfile number.

### 10.4 Multifile: evenly spaced Z with shared X array (XYY)

Similar to 10.3 but with unevenly spaced X values shared across all subfiles.

- Flags: `TMULTI` set, `TXVALS` set, `TXYXYS` clear.
- `TRANDM` and `TORDRD` both clear.
- Interpretation:
  - All subfiles have the same unevenly spaced X values.
  - X values stored once, Y values stored for each subfile.
  - Z spacing determined same as 10.3.

Layout:

1. SPCHDR.
2. Global X array: `fnpts` floats.
3. For each subfile `k` in 0..`fnsub-1`:
   - SUBHDR.
   - Y array: `fnpts` samples.

This format is typically used for multi-spectral kinetics experiments using Raman/CCD/Diode detectors.

### 10.5 Multifile: ordered or random Z

When `TMULTI` is set and:

- `TRANDM` set: arbitrary (random) Z values.
- `TORDRD` set: ordered but uneven Z values.

Then:

- Each subfile's `subtime` is meaningful and must be specified.
- `subnext` is typically equal to `subtime`.
- `fzinc` is ignored; Z values come from each SUBHDR.

For `TORDRD`: Z values must be in ascending or descending order.
For `TRANDM`: Z values can be in any order (though this is not recommended).

The layout is otherwise the same as for evenly spaced Z.

### 10.6 Multifile with per-subfile XY (XYXY)

This is a more complex, but important, mode for data like Mass Spec where each subfile may have a different X grid.

- Flags: `TMULTI` set, `TXVALS` set, `TXYXYS` set.
- Interpretation:
  - Each subfile has its own X array and Y array length given by `subnpts`.
  - An optional SSFSTC directory at `fnpts` provides random access.
  - This is the ONLY format where subfiles can have different point counts.

Layout:

1. SPCHDR.
2. For each subfile:
   - SUBHDR (including `subnpts`).
   - X array: `subnpts` float32s.
   - Y array: `subnpts` samples.
3. Optional SSFSTC directory at byte offset `fnpts` (if `fnpts != 0`).

**Important notes:**

- When `TXYXYS` is set, `fnpts` in the main header is interpreted as:
  - The byte offset to the directory (if non-zero), OR
  - Zero if no directory exists (scan sequentially).
- The X values for each subfile must be in the same orientation (all increasing or all decreasing).
- When determining overall X limits, you must scan through all subfiles and handle both orientations.

This format is typically used for GC-MS data and causes GRAMS to draw "sticks" instead of connected lines.

---

## 11. 4D data (W-plane volumes)

SPC supports 4D data interpreted as volumes in (X, Z, W, Y), where W is an extra axis (e.g. temperature, position, wavelength). This is indicated by:

- `fwplanes` (DWORD) in SPCHDR: number of W-planes.
- `fwinc` (float): W increment when W is evenly spaced.
- `fwtype` (BYTE): W axis unit type.
- `subwlevel` (float) in each SUBHDR: W coordinate for that subfile.

Behavior:

- When `fwplanes != 0`:
  - Subfiles are grouped into planes along W.
  - `fwplanes` gives the number of planes; it must divide the total subfile count `fnsub` evenly.
  - Each plane contains `fnsub / fwplanes` subfiles.
  - Subfiles are "bundled" consecutively: the first N subfiles belong to plane 0, the next N to plane 1, etc.

**Z and W spacing:**

- If `fwinc != 0`, W values are evenly spaced:
  - Starting from `subwlevel` of the first subfile
  - Incremented by `fwinc` after each plane group
  - Formula: $W_p = W_0 + (p \times \text{fwinc})$ where $p$ is plane number

- If `fwinc == 0`, W values are non-evenly spaced:
  - Taken from `subwlevel` of the first subfile in each plane
  - Must be monotonically ordered
  - All subfiles in a plane should have the same `subwlevel`

**Example:**

If `fwplanes = 5` and `fnsub = 20`:
- Each plane has 4 subfiles (20 / 5 = 4)
- Subfiles 0-3 are in W-plane 0
- Subfiles 4-7 are in W-plane 1
- Subfiles 8-11 are in W-plane 2
- etc.
- All subfiles in each group have the same W value

For a minimal reader, you can treat W as a third index over subfiles, using `(plane, z_index_within_plane)`.

**Recommendation**: Equally-spaced W planes (`fwinc != 0`) are recommended as some software may not handle `fwinc == 0`.

---

## 12. Instrument type (`fexper`)

The `fexper` byte describes the data type / instrument technique. In older software, the `TCGRAM` flag in `ftflgs` must be set if `fexper` is non-zero. Common codes:

| Code | Constant    | Meaning                                            |
| ---- | ----------- | -------------------------------------------------- |
| 0    | `SPCGEN`    | General SPC (could be anything)                    |
| 1    | `SPCGC`     | Gas chromatogram                                   |
| 2    | `SPCCGM`    | General chromatogram                               |
| 3    | `SPCHPLC`   | HPLC chromatogram                                  |
| 4    | `SPCFTIR`   | FT-IR / FT-NIR / FT-Raman spectrum or interferogram |
| 5    | `SPCNIR`    | NIR spectra (often calibration datasets)           |
| 7    | `SPCUV`     | UV/VIS spectrum                                    |
| 8    | `SPCXRY`    | X-ray diffraction spectrum                         |
| 9    | `SPCMS`     | Mass spectrum (GC-MS, TOF, etc.)                   |
| 10   | `SPCNMR`    | NMR spectrum or FID                                |
| 11   | `SPCRMN`    | Raman spectrum (not FT-Raman)                      |
| 12   | `SPCFLR`    | Fluorescence spectrum                              |
| 13   | `SPCATM`    | Atomic spectrum                                    |
| 14   | `SPCDAD`    | Chromatography diode array spectra                 |

**Note**: A general chromatogram (`fexper = 0` with `TCGRAM` set) is the same as `SPCGEN`.

This value allows future GRAMS products to be "smart" and tailor operations based on the data type. Many readers simply expose this as metadata and do not special-case per type, except to choose more appropriate default axis labels.

---

## 13. Old 0x4D format (brief)

The original SPC format uses `fversn == 0x4D` and a shorter header type (OSPCHDR). Differences include:

- **Header size**: 256 bytes (not 512)
- `fnpts` is a 32-bit **float** instead of an integer.
- `ffirst` and `flast` are 32-bit **floats** instead of doubles.
- Some fields (`fnsub`, `fmethod`, extra reserved bytes) are missing or different.
- `fexp` is a 16-bit **word** instead of a byte.
- Date/time are stored differently.
- 32-bit Y values have their 16-bit words **swapped** relative to Intel order.
  - Within each word, LSB comes first, but MSW is first overall.
- The first subfile header is included in the main header structure (32 bytes).
- **No support for**: multifiles, Audit Log, XY data, dynamic array sizes.

If you do not need to support very old files (pre-GRAMS/386), you can initially reject files with `fversn != 0x4B` and add 0x4D handling later.

**Recommendation**: All future conversion routines should only create files in the new format (0x4B).

---

## 14. Modification flags (`fmods`)

The `fmods` field in the main header uses bit-encoded flags to track operations applied to the data. Each bit corresponds to a letter representing an operation type:

| Bit | Flag | Operation Type |
|-----|------|----------------|
| 2⁰¹ | A | Averaging (from multiple source traces) |
| 2⁰² | B | Baseline correction or offset functions |
| 2⁰³ | C | Computation (interferogram to spectrum) |
| 2⁰⁴ | D | Derivative (or integrate) functions |
| 2⁰⁶ | E | Resolution Enhancement (deconvolution) |
| 2⁰⁹ | I | Interpolation functions |
| 2¹⁴ | N | Noise reduction smoothing |
| 2¹⁵ | O | Other functions (add, subtract, noise, etc.) |
| 2¹⁹ | S | Spectral Subtraction |
| 2²⁰ | T | Truncation (only portion of original X remains) |
| 2²³ | W | When collected (date/time modified) |
| 2²⁴ | X | X units conversions or X shifting |
| 2²⁵ | Y | Y units conversions (transmission→absorbance, etc.) |
| 2²⁶ | Z | Zap functions (features removed or modified) |

These were widely used in older GRAMS products but have been largely supplanted by the audit trail in the Log Text block. Modern converters typically set this to 0.

---

## 15. Suggested parsing strategy (for a Python/Numpy reader)

A practical approach for an `SPCFile` class:

1. **Read and decode SPCHDR** (512 bytes):
   - Check `fversn` and endianness; currently accept only 0x4B.
   - Extract key fields: `ftflgs`, `fnpts`, `fnsub`, `ffirst`, `flast`, `fexp`, `flogoff`.
   - Parse `fdate` to extract date/time components.

2. **Determine layout** from flags:
   - Single vs multifile: `TMULTI` bit.
   - Implicit vs explicit X: `TXVALS`.
   - Per-subfile X and length: `TXYXYS`.
   - Y data precision: `TSPREC` (16-bit if set).
   - Custom axis labels: `TALABS`.

3. **Read X data**:
   - If no `TXVALS`: compute X from `ffirst`, `flast`, `fnpts`.
   - If `TXVALS` and not `TXYXYS`: read one global X array of `fnpts` floats after SPCHDR.
   - If `TXVALS` and `TXYXYS`: for each subfile, read X arrays of `subnpts` floats after SUBHDR.

4. **Read subfiles sequentially or via directory**:
   - If not `TXYXYS`: iterate subfiles in order, each with SUBHDR followed by Y data.
   - If `TXYXYS`:
     - If `fnpts != 0`: optionally read SSFSTC directory for random access.
     - If `fnpts == 0`: scan sequentially through subfiles.

5. **Decode Y arrays**:
   - Determine per-subfile exponent (`subexp` if `TMULTI`, else `fexp`).
   - Check if exponent is 0x80 (float Y data).
   - Determine Y sample width from `TSPREC`.
   - Read raw integers/floats and scale to Python floats (or Numpy arrays).
   - For fixed-point: apply formula $y = 2^{exp} \cdot ival / 2^N$

6. **Parse Z axis**:
   - For evenly-spaced multifiles: use `fzinc` or first subfile's `subnext - subtime`.
   - For ordered/random: read `subtime` from each SUBHDR.

7. **Handle 4D data** (if `fwplanes != 0`):
   - Group subfiles into W-planes.
   - Extract W values from `subwlevel` and/or `fwinc`.

8. **Optional: parse log block** if `flogoff != 0`:
   - Read LOGSTC header.
   - Skip binary blocks (`logbins`, `logdsks`).
   - Extract log text at offset `logtxto`.
   - For NMR files: read imaginary data from binary block.

This strategy keeps the implementation manageable while covering single-spectrum, multifile, and XY / TXYXYS cases.

---

## 16. GSPCIO library (reference only)

The GSPCIO library is a COM-based component provided by Thermo Galactic for reading/writing SPC files in Windows applications. While the Brief Guide PDF focuses heavily on this library, it is not necessary for implementing a file parser. Key points:

- **Language**: Designed for Visual Basic / COM applications.
- **Functionality**: Provides high-level methods to read/write SPC files without dealing with binary structures directly.
- **Modern use**: Most modern parsers implement direct binary reading rather than using GSPCIO.

GSPCIO key methods (for reference only):

- `CreateSPC`: Create new SPC file in memory
- `OpenFile`: Open SPC file from disk
- `OpenBlob`: Open SPC file from memory
- `SaveFile`: Save SPC file to disk
- `AddPoints`: Add Y (and optionally X) data
- `AddSubfile`: Add a subfile to a multifile
- Property accessors: `FirstX`, `LastX`, `NumPoints`, `XPoints`, `YPoints`, etc.

**Note**: GSPCIO always writes Y data in 32-bit IEEE floating point format, even though it can read all formats. This is the recommended approach for new converters.

---

## 17. Out-of-scope details

The SPC SDK and UDF specification include additional topics not fully covered here:

- Detailed GRAMS DDE structures and file converter APIs (`gcdll.h`, `gcdllerr.h`).
- Full semantics of post-processing codes (`fprocs`) and related workflow metadata.
- Technique-specific log text parameters (FT-IR, NIR, UV/VIS, NMR, MS, chromatography).
- GLP compliance and audit trail requirements.
- Calibration File List (CFL) format for NIR quantitative analysis.
- Complete list of modification flags and their historical usage.
- Legacy behavior of the old 0x4D header beyond what is summarized.

For those, refer directly to the original headers and documentation as needed.

---

## 18. Best practices for file writers

When creating new SPC files, follow these recommendations:

1. **Use the new format**: Always write `fversn = 0x4B` (little-endian) format.

2. **Y data format**: Write Y data as 32-bit IEEE floats (`fexp = 0x80`) unless higher precision is required by the source data format.

3. **Set all required fields**:
   - Fill in `ffirst`, `flast` even for XY files.
   - Set `subindx` correctly for all subfiles (0-indexed).
   - Set `fdate` to the actual collection date/time.

4. **Axis labels**: Use standard enumerated types where possible rather than custom labels.

5. **Log text**: Include common parameters like `MODEL`, `MEMO`, `USER`, acquisition parameters.

6. **Z spacing**: For evenly-spaced multifiles, set `fzinc` explicitly rather than relying on the first subfile calculation.

7. **W-planes**: If using 4D data, prefer evenly-spaced W-planes (`fwinc != 0`).

8. **TXYXYS directory**: When creating XYXY multifiles, always include the directory unless the file will only be accessed sequentially.

9. **Reserved fields**: Always zero out reserved fields and spare areas.

10. **File extensions**: Use appropriate extensions (`.SPC` for spectra, `.CGM` for chromatograms).

---

## Appendix A: Complete structure sizes

| Structure | Size (bytes) | Purpose |
|-----------|-------------|---------|
| SPCHDR    | 512         | Main file header (new format) |
| OSPCHDR   | 256         | Main file header (old 0x4D format) |
| SUBHDR    | 32          | Subfile header |
| SSFSTC    | 12          | XY directory entry |
| LOGSTC    | 64          | Log block header |

---

## Appendix B: Useful constants

```python
# ftflgs bits
TSPREC = 0x01  # 16-bit Y data
TCGRAM = 0x02  # Chromatogram / fexper enabled
TMULTI = 0x04  # Multifile
TRANDM = 0x08  # Random Z values
TORDRD = 0x10  # Ordered but uneven Z
TALABS = 0x20  # Custom axis labels
TXYXYS = 0x40  # Per-subfile X arrays
TXVALS = 0x80  # X values stored

# fversn values
FVERSN_NEW_LSB = 0x4B  # New format, little-endian
FVERSN_NEW_MSB = 0x4C  # New format, big-endian
FVERSN_OLD = 0x4D      # Old format

# Y data float flag
FLOAT_EXP = 0x80  # -128 signed

# subflgs bits
SUBCHGD = 0x01  # Subfile changed
SUBNOPT = 0x08  # Don't use peak table
SUBMODF = 0x80  # Modified by arithmetic
```

---

## Appendix C: Common pitfalls

1. **Forgetting SUBHDR**: Even single-spectrum files have a SUBHDR before Y data.

2. **Wrong exponent**: Using `fexp` for multifiles instead of per-subfile `subexp`.

3. **TXYXYS fnpts**: Forgetting that `fnpts` is a byte offset, not a point count, in TXYXYS mode.

4. **16-bit Y data**: Not checking `TSPREC` flag before reading Y values.

5. **Z increment**: Not handling the case where `fzinc = 0` and you must use the first subfile's `subnext - subtime`.

6. **Float Y detection**: Not checking for `exp == 0x80` before applying fixed-point formula.

7. **X array position**: Reading X array in the wrong location (after SPCHDR for global, after SUBHDR for per-subfile).

8. **Log text offset**: Forgetting that `logtxto` is relative to the start of LOGSTC, not the file.

9. **Custom labels**: Not checking `TALABS` flag before using `fcatxt`.

10. **4D grouping**: Not dividing subfiles correctly into W-planes.

---

## Appendix D: Quick decision tree

```
┌─ fversn == 0x4B? ─┐
│  No: Old format   │
│  Yes: Continue    │
└───────────────────┘
       │
┌─ TMULTI set? ─────┐
│  No: Single file  │
│  Yes: Multifile   │
└───────────────────┘
       │
┌─ TXVALS set? ─────┐
│  No: Y-only       │
│  Yes: XY data     │
└───────────────────┘
       │
┌─ TXYXYS set? ─────┐
│  No: Shared X     │
│  Yes: Per-sub 