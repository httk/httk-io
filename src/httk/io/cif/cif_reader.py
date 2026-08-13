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

import logging
import os
import re
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any, Self

from httk.core import TextstreamFileView

from .cif_tags import CIF_TAGS

logger = logging.getLogger(__name__)


# Normalized tags consumed by cif_parser.py and mcif_parser.py; a loop is droppable only
# when none of its columns is in this set. The Fourier coefficient template is checked below.
_PROTECTED_LOOP_TAGS = frozenset(
    {
        'cell_modulation_dimension',
        'cell_length_a',
        'cell_length_b',
        'cell_length_c',
        'cell_angle_alpha',
        'cell_angle_beta',
        'cell_angle_gamma',
        'atom_site_type_symbol',
        'atom_site_label',
        'atom_site_fract_x',
        'atom_site_fract_y',
        'atom_site_fract_z',
        'atom_site_occupancy',
        'atom_site_wyckoff_label',
        'atom_site_symmetry_multiplicity',
        'httk_atom_site_fract_x_exact',
        'httk_atom_site_fract_y_exact',
        'httk_atom_site_fract_z_exact',
        'httk_atom_site_occupancy_exact',
        'space_group_symop.operation_xyz',
        'space_group_symop_operation_xyz',
        'symmetry_equiv_pos_as_xyz',
        'space_group_name_h-m_alt',
        'symmetry_space_group_name_h-m',
        'space_group_name_hall',
        'symmetry_space_group_name_hall',
        'space_group_it_number',
        'symmetry_space_group_it_number',
        'symmetry_int_tables_number',
        'database_code_icsd',
        'citation_doi',
        'parent_propagation_vector.kxkykz',
        'atom_site_moment.label',
        'atom_site_moment.crystalaxis_x',
        'atom_site_moment.crystalaxis_y',
        'atom_site_moment.crystalaxis_z',
        'space_group_symop_magn_operation.xyz',
        'space_group_symop_magn_ssg_operation.algebraic',
        'space_group_symop_magn_centering.xyz',
        'space_group_symop_magn_ssg_centering.algebraic',
        'space_group_magn.number_bns',
        'space_group_magn.name_bns',
        'parent_space_group.name_h-m_alt',
        'parent_space_group.it_number',
        *CIF_TAGS['structural_q'],
        CIF_TAGS['structural_displacement_label'],
        CIF_TAGS['structural_occupancy_label'],
        *CIF_TAGS['magnetic_cartesian_moment'],
        CIF_TAGS['magnetic_fourier_label'],
        CIF_TAGS['magnetic_ssg_name'],
    }
)
_MAGNETIC_FOURIER_COEFFICIENT_RE = re.compile(
    '^' + re.escape(CIF_TAGS['magnetic_fourier_coeff']).replace(r'\{\}', r'\d+') + '$'
)
_STRUCTURAL_LOOP_PREFIXES = ('atom_site', 'space_group', 'symmetry')


def _is_consumed_tag(name: str) -> bool:
    return (
        name in _PROTECTED_LOOP_TAGS
        or name.startswith(_STRUCTURAL_LOOP_PREFIXES)
        or _MAGNETIC_FOURIER_COEFFICIENT_RE.fullmatch(name) is not None
    )


def _is_repairable_loop(header: list[str]) -> bool:
    return not any(_is_consumed_tag(name) for name in header)


class _RewindableIterator:
    def __init__(self, iterator: Iterable[str]) -> None:
        self._iter: Iterator[str] = iter(iterator)
        self._rewind = False
        self._cache: str | None = None

    def __iter__(self) -> Self:
        return self

    def __next__(self) -> str:
        if self._rewind:
            self._rewind = False
        else:
            self._cache = next(self._iter)
        assert self._cache is not None
        return self._cache

    def rewind(self, rewindstr: str | None = None) -> None:
        if self._rewind:
            raise RuntimeError("Tried to backup more than one step.")
        elif self._cache is None:
            raise RuntimeError("Can't backup past the beginning.")
        self._rewind = True
        if rewindstr is not None:
            self._cache = rewindstr


def _read_cif_rewind_if_needed(f: _RewindableIterator, row: str, done_fields: int) -> bool:
    splitstr = row.lstrip().split(None, done_fields)
    if len(splitstr) > 1:
        rest = splitstr[-1]
        if rest.strip() != "":
            f.rewind(rest)
            return True
        return False
    else:
        return False


def _read_cif_loop(
    f: _RewindableIterator,
    pragmatic: bool = True,
    allow_cif2: bool = False,
    *,
    block_name: str,
    autocorrect: bool = False,
) -> dict[str, list[Any]] | None:
    noteol = False
    loop_data: dict[str, list[Any]] = {}
    header = []
    for row in f:
        striprow = row.strip()
        lowrow = striprow.lower()
        if lowrow.startswith("_"):
            loop_data[lowrow[1:]] = []
            header += [lowrow[1:]]
            noteol = _read_cif_rewind_if_needed(f, row, 1)
        else:
            f.rewind()
            break

    while len(header) > 0:
        for i in range(len(header)):
            try:
                row = next(f)
                while row.isspace():
                    row = next(f)
            except StopIteration:
                break
            striprow = row.strip()
            lowrow = striprow.lower()
            if not row or row.startswith("_") or lowrow.startswith(("data_", "loop_")):
                f.rewind()
                break
            f.rewind()
            val, noteol = _read_cif_data_value(f, noteol, pragmatic, allow_cif2, inloop=True)
            if val is None:
                # Could be a comment line, etc.
                continue
            loop_data[header[i]].append(val)
        else:
            continue
        break
    counts = {name: len(values) for name, values in loop_data.items()}
    if len(set(counts.values())) > 1:
        rendered_counts = ", ".join(f"{name}={count}" for name, count in counts.items())
        message = f"CIF loop with {len(header)} columns has mismatched value counts: {rendered_counts}"
        if _is_repairable_loop(header):
            if autocorrect:
                logger.warning(
                    "CIF block %r: dropped malformed auxiliary loop starting with _%s",
                    block_name,
                    header[0],
                    extra={'context': 'cif'},
                )
                return None
            message += " (an auxiliary loop like this can be dropped by loading with autocorrect=True, which applies documented repairs with warnings)"
        raise ValueError(message)
    return loop_data


def _read_cif_data_value(
    f: _RewindableIterator,
    noteol: bool,
    pragmatic: bool = True,
    allow_cif2: bool = False,
    inloop: bool = False,
    inlist: bool = False,
) -> tuple[Any, bool]:
    data_value: Any = None
    for row in f:
        striprow = row.strip()
        if striprow.startswith("#") or striprow == "":
            noteol = False
            continue
        elif (not noteol) and row.startswith(';'):
            folded = False
            newline = False
            data_value = ""
            if row[1] == "\\" and row[2:].rstrip("\r\n") == "":
                folded = True
            elif row[1:].isspace():
                if not pragmatic:
                    data_value = row.lstrip().rstrip('\r\n')
                    newline = True
            else:
                data_value = row.lstrip()[1:].rstrip('\r\n')
                newline = True
            stripirow = ""
            for irow in f:
                stripirow = irow.strip()
                if irow.startswith(';'):
                    break
                if newline:
                    data_value += '\n'
                    newline = True
                if folded and irow.rstrip('\r\n').endswith("\\"):
                    data_value += irow.rstrip('\r\n').rstrip("\\")
                    newline = False
                else:
                    data_value += irow.rstrip('\r\n')
                    newline = True
            if len(stripirow) > 1:
                f.rewind(stripirow[1:])
                noteol = True
            else:
                noteol = False
            break
        elif striprow.startswith(("'", '"')):
            # The cif quoting rules are ... weird. Quotes are "escaped" if they are not followed by whitespace.
            quote = striprow[0]
            starti = 1
            for chari in range(1, len(striprow) - 1):
                if striprow[chari] == quote and str(striprow[chari + 1]).isspace():
                    endi = chari
                    endq = chari + 1
                    break
            else:
                if striprow[-1] != quote:
                    starti = 0
                    endi = len(striprow)
                    endq = len(striprow)
                else:
                    endi = len(striprow) - 1
                    endq = len(striprow)
            data_value = striprow[starti:endi]
            if endq != len(striprow):
                f.rewind(striprow[endq:])
                noteol = True
            else:
                noteol = False
            break
        elif allow_cif2 and inlist and striprow.startswith("]"):
            # TODO: Is ] allowed without whitespace after? I need to check the spec
            splitstr = striprow.split("]", 1)
            if len(splitstr) > 1 and len(splitstr[1]) > 0:
                f.rewind(splitstr[1])
                noteol = True
            data_value = None
            break
        elif allow_cif2 and striprow.startswith("["):
            if len(striprow) > 1:
                f.rewind(striprow[1:])
                noteol = True
            data_value = []
            while True:
                innerval, noteol = _read_cif_data_value(f, noteol, pragmatic, allow_cif2, inloop=False, inlist=True)
                if innerval is None:
                    break
                data_value += [innerval]
            break
        elif allow_cif2 and inlist and ("]" in striprow):
            splitstr2 = striprow.split("]", 1)
            splitstr = splitstr2[0].split(None, 1)
            data_value = splitstr[0].strip()
            rightside = ""
            if len(splitstr) > 1:
                f.rewind(splitstr[1] + "]" + splitstr2[1])
            else:
                f.rewind("]" + splitstr2[1])
            noteol = True
            break
        else:
            if pragmatic and not inloop:
                # In pragmatic mode, if we are not in a loop and there is more than one data value
                # separated by whitespace, read all of it. This should always be ok to do, since
                # multiple data values in this situation would be an
                # error in the file otherwise, but if there is whitespace + underscore/data_/loop_ we parse that
                # as a new symbol, since otherwise we COULD misread valid files (with very weird formatting...).
                splitstr = re.split(r'\s+_|\s+data_|\s+loop_', striprow, maxsplit=1)
            else:
                splitstr = striprow.split(None, 1)
            # "Data on a line following a hash character `#' is considered to be a comment,
            # except if it is contained within a text string."
            data_value = splitstr[0].partition("#")[0].strip()
            rightside = ""
            if len(splitstr) > 1:
                rightside = splitstr[1].strip()
            if rightside != "":
                f.rewind(rightside)
                noteol = True
            else:
                noteol = False
            break
    return data_value, noteol


def _read_cif_data_block(
    f: _RewindableIterator,
    pragmatic: bool = True,
    allow_cif2: bool = False,
    *,
    block_name: str,
    autocorrect: bool = False,
) -> dict[str, Any]:
    data_items: dict[str, Any] = {}
    loops = 0
    for row in f:
        striprow = row.strip()
        lowrow = striprow.lower()
        if striprow.startswith("#"):
            continue
        elif lowrow.startswith("data_"):
            f.rewind()
            return data_items
        elif lowrow.startswith("loop_"):
            _read_cif_rewind_if_needed(f, row, 1)
            loopdata = _read_cif_loop(f, pragmatic, allow_cif2, block_name=block_name, autocorrect=autocorrect)
            if loopdata is None:
                continue
            data_items['loop_' + str(loops)] = list(loopdata.keys())
            loops += 1
            data_items.update(loopdata)
        elif striprow.startswith(";"):
            # Multi-line string that we've failed to tie to a name, lets just skip it, maybe we should warn
            for irow in f:
                if irow.rstrip() == ";":
                    break
        elif striprow.startswith("_"):
            lowsplit = lowrow.split()
            data_name = lowsplit[0][1:]
            if len(lowsplit) > 1:
                noteol = True
                rightside = striprow.split(None, 1)[1].strip()
                f.rewind(rightside)
            else:
                noteol = False
            data_value, noteol = _read_cif_data_value(f, noteol, pragmatic, allow_cif2, inloop=False)
            data_items[data_name] = data_value
    return data_items


def _read_cif(
    f: _RewindableIterator, pragmatic: bool, allow_cif2: bool, autocorrect: bool
) -> tuple[list[tuple[str, dict[str, Any]]], str]:
    header = ""
    datalist = []
    for row in f:
        if row.strip().startswith("#"):
            header += row
        else:
            f.rewind()
            break

    for row in f:
        lowrow = row.strip().lower()
        if lowrow.startswith("data_"):
            data_block_name = lowrow.partition('_')[2].split()[0].strip()
            _read_cif_rewind_if_needed(f, row, 1)
            datalist.append(
                (
                    data_block_name,
                    _read_cif_data_block(f, pragmatic, allow_cif2, block_name=data_block_name, autocorrect=autocorrect),
                )
            )
    return datalist, header


def read_cif(
    source: str | os.PathLike[str] | Iterable[str],
    pragmatic: bool = True,
    allow_cif2: bool = False,
    *,
    autocorrect: bool = False,
) -> tuple[list[tuple[str, dict[str, Any]]], str]:
    """Read CIF text as ``(data_blocks, header)``.

    Paths are opened through :class:`httk.core.TextstreamFileView`, including compressed
    CIF files. Open streams and iterables are consumed but left open.

    :param source: A filename, open text stream, or iterable of CIF lines.
    :param pragmatic: Accept selected common deviations from strict CIF tokenization.
    :param allow_cif2: Parse CIF2 list values in addition to CIF1 data.
    :param autocorrect: Drop malformed auxiliary loops and warn about each repair.
    :return: The data blocks and the leading comment header.
    :raises ValueError: If a loop contains mismatched column value counts.
    """
    if isinstance(source, (str, os.PathLike)):
        with TextstreamFileView(Path(source)) as stream:
            return _read_cif(_RewindableIterator(stream), pragmatic, allow_cif2, autocorrect)
    return _read_cif(_RewindableIterator(source), pragmatic, allow_cif2, autocorrect)
