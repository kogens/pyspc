from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .parser import read_all_subheaders, read_header


@dataclass
class SPCSubfile:
    """Single subfile (spectrum) in an SPC file."""

    x: np.ndarray
    y: np.ndarray


class SPCFile:
    """In-memory representation of a GRAMS SPC file."""

    def __init__(self, path: str | Path) -> None:
        self._path = Path(path)

        if not self._path.is_file():
            raise FileNotFoundError(self._path)
        
        with self._path.open("rb") as f:
            self.header = read_header(f)

        if self.header['version'] != 0x4B:
            raise ValueError(f"Unsupported SPC version: {self.header['version']:02X}")
        
        # Read all subheaders using parser
        with self._path.open("rb") as f:
            self._raw_subheaders = read_all_subheaders(f, self.header)
        
        # Temporary stub data so tests can exercise the public API
        self._subfiles: list[SPCSubfile] = [
            SPCSubfile(x=np.array([], dtype=float), y=np.array([], dtype=float))
        ]

    @property
    def path(self) -> Path:
        """Filesystem path of the underlying SPC file."""

        return self._path

    @property
    def subfiles(self) -> list[SPCSubfile]:
        """All subfiles contained in this SPC file."""

        return self._subfiles

    @property
    def x(self) -> np.ndarray:
        """X values for the file.

        For single-spectrum files this is the X axis of the only
        subfile. This behavior will be refined as multifile support
        is implemented.
        """

        return self._subfiles[0].x

    @property
    def y(self) -> np.ndarray:
        """Y values for the file.

        For single-spectrum files this is the Y axis of the only
        subfile. This behavior will be refined as multifile support
        is implemented.
        """

        return self._subfiles[0].y

    def __len__(self) -> int:
        return len(self._subfiles)

    def __getitem__(self, index: int) -> SPCSubfile:
        return self._subfiles[index]

    def __iter__(self) -> Iterable[SPCSubfile]:
        return iter(self._subfiles)
