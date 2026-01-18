"""Data structures for SPC files."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class SPCSubfile:
    """Single subfile (spectrum) in an SPC file."""

    x: np.ndarray
    y: np.ndarray
    header: dict[str, object] | None = None

    @classmethod
    def from_raw(
        cls, x: np.ndarray, y: np.ndarray, header: dict[str, object] | None = None
    ) -> SPCSubfile:
        """Create an SPCSubfile from raw X/Y arrays and optional subheader.

        Args:
            x: X coordinate array
            y: Y value array
            header: Optional subheader dict from parser

        Returns:
            New SPCSubfile instance
        """
        return cls(x=x, y=y, header=header)
