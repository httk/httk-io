#
#    The high-throughput toolkit (httk)
#    Copyright (C) 2012-2025 The httk AUTHORS
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""A string-preserving reader for VASP POSCAR/CONTCAR files.

:func:`read_poscar` parses a POSCAR/CONTCAR file into a neutral, JSON-able
mapping whose numeric fields are kept as the **verbatim strings** found in the
file. It performs no numeric conversion and imports nothing from
*httk-atomistic*; turning the mapping into a ``Structure`` is the job of
``httk.atomistic.structure_from_poscar``.
"""

import os
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from httk.core import TextstreamFileView


def _next_line(lines: Iterator[str], lineno: int, what: str) -> str:
    """Return the next line's text, raising a clear error at end of file."""
    try:
        return next(lines)
    except StopIteration:
        raise ValueError(f"Malformed POSCAR: unexpected end of file at line {lineno} (expected {what}).") from None


def _parse_flag(token: str, lineno: int) -> bool:
    if token in ("T", "t"):
        return True
    if token in ("F", "f"):
        return False
    raise ValueError(f"Malformed POSCAR line {lineno}: expected a selective-dynamics flag 'T' or 'F', got {token!r}.")


def read_poscar(source: Any) -> dict[str, Any]:
    """Parse a VASP POSCAR/CONTCAR into a neutral, string-preserving mapping.

    ``source`` may be a filename (``str`` or :class:`os.PathLike`, opened through
    :class:`httk.core.TextstreamFileView` so compressed files such as
    ``CONTCAR.bz2`` are decompressed transparently) or an already-open text
    stream / iterable of lines.

    The returned mapping has the keys ``format`` (always ``"vasp-poscar"``),
    ``comment``, ``scale`` and ``volume`` (exactly one is a string, the other
    ``None``), ``cell`` (a 3x3 list of coordinate strings), ``symbols``
    (``list[str]`` for VASP-5, ``None`` for VASP-4), ``counts`` (``list[int]``),
    ``cartesian`` (``bool``), ``coords`` (an ``N``x3 list of strings), and
    ``selective_dynamics`` (an ``N``x3 list of booleans when the file declares
    selective dynamics, else ``None``). Malformed input raises a clear
    :class:`ValueError` naming the offending line.
    """
    opened: TextstreamFileView | None = None
    if isinstance(source, (str, os.PathLike)):
        opened = TextstreamFileView(Path(source))
        raw: Iterable[str] = opened
    else:
        raw = source
    try:
        return _read_poscar(iter(raw))
    finally:
        if opened is not None:
            opened.close()


def _read_poscar(lines: Iterator[str]) -> dict[str, Any]:
    comment = _next_line(lines, 1, "comment").strip()

    scale_line = _next_line(lines, 2, "scale/volume").strip()
    try:
        scale_value = float(scale_line)
    except ValueError:
        raise ValueError(f"Malformed POSCAR line 2: scale/volume {scale_line!r} is not a number.") from None
    if scale_value < 0:
        # A negative universal scaling factor means |value| is the target VOLUME.
        volume: str | None = scale_line[1:] if scale_line.startswith("-") else scale_line
        scale: str | None = None
    else:
        scale = scale_line
        volume = None

    cell: list[list[str]] = []
    for i in range(3):
        lineno = 3 + i
        tokens = _next_line(lines, lineno, "a lattice-vector row").strip().split()
        if len(tokens) < 3:
            raise ValueError(f"Malformed POSCAR line {lineno}: expected 3 lattice-vector components, got {tokens!r}.")
        cell.append(tokens[:3])

    species_line = _next_line(lines, 6, "species symbols or atom counts").strip().split()
    if not species_line:
        raise ValueError("Malformed POSCAR line 6: expected species symbols or atom counts, got a blank line.")
    try:
        counts = [int(token) for token in species_line]
        symbols: list[str] | None = None
        counts_lineno = 6
    except ValueError:
        symbols = species_line
        counts_line = _next_line(lines, 7, "atom counts").strip().split()
        try:
            counts = [int(token) for token in counts_line]
        except ValueError:
            raise ValueError(f"Malformed POSCAR line 7: atom counts {counts_line!r} are not all integers.") from None
        if len(counts) != len(symbols):
            raise ValueError(f"Malformed POSCAR: {len(symbols)} species symbol(s) but {len(counts)} atom count(s).")
        counts_lineno = 7

    n_atoms = sum(counts)

    # Optional selective-dynamics line, then the coordinate-type line.
    mode_lineno = counts_lineno + 1
    mode_line = _next_line(lines, mode_lineno, "coordinate type (or 'Selective dynamics')").strip()
    selective = bool(mode_line) and mode_line[0] in "Ss"
    if selective:
        coordtype_lineno = mode_lineno + 1
        coordtype = _next_line(lines, coordtype_lineno, "coordinate type").strip()
    else:
        coordtype_lineno = mode_lineno
        coordtype = mode_line
    if not coordtype:
        raise ValueError(f"Malformed POSCAR line {coordtype_lineno}: missing coordinate type (Direct/Cartesian).")
    cartesian = coordtype[0] in "CcKk"

    coords: list[list[str]] = []
    selective_dynamics: list[list[bool]] | None = [] if selective else None
    for i in range(n_atoms):
        lineno = coordtype_lineno + 1 + i
        tokens = _next_line(lines, lineno, "an atomic coordinate row").strip().split()
        if len(tokens) < 3:
            raise ValueError(f"Malformed POSCAR line {lineno}: expected 3 coordinate components, got {tokens!r}.")
        coords.append(tokens[:3])
        if selective_dynamics is not None:
            if len(tokens) < 6:
                raise ValueError(
                    f"Malformed POSCAR line {lineno}: selective dynamics declared but flags are missing ({tokens!r})."
                )
            selective_dynamics.append([_parse_flag(tokens[3 + j], lineno) for j in range(3)])

    return {
        "format": "vasp-poscar",
        "comment": comment,
        "scale": scale,
        "volume": volume,
        "cell": cell,
        "symbols": symbols,
        "counts": counts,
        "cartesian": cartesian,
        "coords": coords,
        "selective_dynamics": selective_dynamics,
    }
