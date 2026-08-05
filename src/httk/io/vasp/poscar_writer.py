"""Writer for the neutral, string-preserving VASP POSCAR payload."""

import io
import os
from collections.abc import Mapping
from contextlib import ExitStack
from typing import Any


def _write_poscar_payload(destination: str | os.PathLike[str] | io.TextIOBase, data: Mapping[str, Any]) -> None:
    """Write a neutral payload returned by :func:`read_poscar`."""
    if "format" in data and data["format"] != "vasp-poscar":
        raise ValueError("POSCAR payload format must be 'vasp-poscar'.")

    scale = data.get("scale")
    volume = data.get("volume")
    if isinstance(scale, str) and volume is None:
        scale_line = scale
    elif scale is None and isinstance(volume, str):
        scale_line = "-" + volume
    else:
        raise ValueError("POSCAR payload must have exactly one of scale or volume as a non-None string.")

    comment = data.get("comment")
    if not isinstance(comment, str):
        raise ValueError("POSCAR payload comment must be a string.")
    cell = data.get("cell")
    if (
        not isinstance(cell, list)
        or len(cell) != 3
        or any(
            not isinstance(row, list) or len(row) != 3 or any(not isinstance(token, str) for token in row)
            for row in cell
        )
    ):
        raise ValueError("POSCAR payload cell must be 3 rows of 3 string tokens.")

    counts = data.get("counts")
    if not isinstance(counts, list) or not counts or any(type(count) is not int or count < 0 for count in counts):
        raise ValueError("POSCAR payload counts must be a non-empty list of non-negative integers.")
    n_atoms = sum(counts)

    coords = data.get("coords")
    if (
        not isinstance(coords, list)
        or len(coords) != n_atoms
        or any(
            not isinstance(row, list) or len(row) != 3 or any(not isinstance(token, str) for token in row)
            for row in coords
        )
    ):
        raise ValueError(f"POSCAR payload coords must have exactly {n_atoms} rows of 3 string tokens.")

    symbols = data.get("symbols")
    if symbols is not None and (
        not isinstance(symbols, list)
        or len(symbols) != len(counts)
        or any(not isinstance(symbol, str) for symbol in symbols)
    ):
        raise ValueError("POSCAR payload symbols must match the number of count entries and contain strings.")

    selective_dynamics = data.get("selective_dynamics")
    if selective_dynamics is not None and (
        not isinstance(selective_dynamics, list)
        or len(selective_dynamics) != n_atoms
        or any(
            not isinstance(row, list) or len(row) != 3 or any(type(flag) is not bool for flag in row)
            for row in selective_dynamics
        )
    ):
        raise ValueError(f"POSCAR payload selective_dynamics must have {n_atoms} rows of 3 booleans.")

    cartesian = data.get("cartesian")
    if not isinstance(cartesian, bool):
        raise ValueError("POSCAR payload cartesian must be a boolean.")

    with ExitStack() as stack:
        f: io.TextIOBase
        if isinstance(destination, (str, os.PathLike)):
            f = stack.enter_context(open(destination, "w", encoding="utf-8"))
        else:
            f = destination

        f.write(comment + "\n")
        f.write(scale_line + "\n")
        for row in cell:
            f.write(" ".join(row) + "\n")
        if symbols is not None:
            f.write(" ".join(symbols) + "\n")
        f.write(" ".join(str(count) for count in counts) + "\n")
        if selective_dynamics is not None:
            f.write("Selective dynamics\n")
        f.write(("Cartesian" if cartesian else "Direct") + "\n")
        for row, flags in zip(coords, selective_dynamics or []):
            f.write(" ".join(row + ["T" if flag else "F" for flag in flags]) + "\n")
        if selective_dynamics is None:
            for row in coords:
                f.write(" ".join(row) + "\n")
