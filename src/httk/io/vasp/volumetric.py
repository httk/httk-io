"""Writer for VASP/VESTA volumetric grids."""

import io
import os
from collections.abc import Mapping
from contextlib import ExitStack
from typing import Any

try:
    import numpy
except ImportError:
    raise ImportError("httk.io.vasp.volumetric requires numpy; install httk-io[numpy]") from None

from .poscar_writer import _write_poscar_payload


def write_vasp_volumetric(
    destination: str | os.PathLike[str] | io.TextIOBase,
    poscar_payload: Mapping[str, Any],
    grid: Any,
    *,
    cols: int = 10,
) -> None:
    """Write POSCAR content followed by a Fortran-order real-valued grid.

    :param destination: Filesystem path or open text stream for the output.
    :param poscar_payload: Neutral POSCAR payload supplying the structure header.
    :param grid: Three-dimensional real-valued grid to write in Fortran order.
    :param cols: Maximum number of values written on each grid line.
    :raises ValueError: If ``cols`` or ``grid`` cannot be represented in the format.
    """
    if isinstance(cols, bool) or not isinstance(cols, int) or cols <= 0:
        raise ValueError("volumetric cols must be a positive integer.")
    grid = numpy.asarray(grid)
    if grid.ndim != 3:
        raise ValueError("volumetric grid must be a three-dimensional array.")
    if numpy.iscomplexobj(grid):
        raise ValueError("volumetric grid must be real; pass grid.real or grid.imag explicitly.")

    with ExitStack() as stack:
        file: io.TextIOBase
        if isinstance(destination, (str, os.PathLike)):
            file = stack.enter_context(open(destination, "w", encoding="utf-8"))
        else:
            file = destination
        _write_poscar_payload(file, poscar_payload)
        file.write("\n")
        file.write("{} {} {}\n".format(*grid.shape))
        values = numpy.asarray(grid, dtype=numpy.float64).flatten(order="F")
        for start in range(0, len(values), cols):
            file.write(" ".join(f"{value:16.8E}" for value in values[start : start + cols]) + "\n")
