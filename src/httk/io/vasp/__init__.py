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

"""VASP file-format readers for *httk-io*."""

import io
import os
from collections.abc import Mapping
from typing import Any

from .oszicar import read_oszicar
from .outcar import ElasticModuliBlock, OutcarFile, OutcarFrame, read_outcar
from .outputs import VASPOutputs
from .poscar_reader import read_poscar
from .potcar import read_potcar_summary
from .xdatcar import XdatcarFile, read_xdatcar


def read_wavecar(
    source: str | os.PathLike[str], *, double_precision: bool | None = None, gamma_half: str | None = None
) -> dict[str, Any]:
    """Lazily read a VASP WAVECAR path into a neutral payload."""
    from .wavecar import read_wavecar as _read_wavecar

    return _read_wavecar(source, double_precision=double_precision, gamma_half=gamma_half)


def write_wavecar(destination: str | os.PathLike[str], payload: Mapping[str, Any]) -> None:
    """Lazily write a neutral WAVECAR payload to a binary path."""
    from .wavecar import write_wavecar as _write_wavecar

    _write_wavecar(destination, payload)


def write_vasp_volumetric(
    destination: str | os.PathLike[str] | io.TextIOBase,
    poscar_payload: Mapping[str, Any],
    grid: Any,
    *,
    cols: int = 10,
) -> None:
    """Lazily write POSCAR content and a VASP/VESTA volumetric grid."""
    from .volumetric import write_vasp_volumetric as _write_vasp_volumetric

    _write_vasp_volumetric(destination, poscar_payload, grid, cols=cols)


__all__ = [
    "ElasticModuliBlock",
    "OutcarFile",
    "OutcarFrame",
    "VASPOutputs",
    "XdatcarFile",
    "read_oszicar",
    "read_outcar",
    "read_poscar",
    "read_potcar_summary",
    "read_wavecar",
    "read_xdatcar",
    "write_vasp_volumetric",
    "write_wavecar",
]
