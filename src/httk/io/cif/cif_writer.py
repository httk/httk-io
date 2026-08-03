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

import io
import os
import re
import sys
import warnings
from collections.abc import Iterable, Mapping, Sequence
from contextlib import ExitStack
from decimal import Decimal, localcontext
from fractions import Fraction
from typing import Any, cast

_cif_ordinary_char = r'!%&()*+,-./0123456789:<=>?@ABCDEFGHIHJKLMNOPQRSTUVWXYZ\^`abcdefghijklmnopqrstuvwxyz{|}~'
_cif_non_blank_char = _cif_ordinary_char + '"' + "#$" + "'" + "_" + ";[]"
_cif_text_lead_char = _cif_ordinary_char + '"' + "#$" + "'" + "_ \t[]"
_cif_any_print_char = _cif_ordinary_char + '"' + "#$" + "'" + "_ \t;[]"
_cif_non_blank_char_table = str.maketrans(_cif_non_blank_char, _cif_non_blank_char)
_cif_unicode_translation_table: dict[int, int | None] = {}
for i in range(sys.maxunicode + 1):
    _cif_unicode_translation_table[i] = None
for key, value in _cif_non_blank_char_table.items():
    _cif_unicode_translation_table[key] = value
_cif_integer_regex = re.compile(r'^[+-]?[0-9]+$')
_cif_float_regex = re.compile(r'^[+-]?[0-9]+[eE][+-]?[0-9]+|([+-]?[0-9]*\.[0-9]+|[+-]?[0-9]\.)([eE][+-]?[0-9]+)?$')
_cif_simplestring_regex = re.compile(r'^[A-Za-z0-9()][A-Za-z0-9()+-]*$')
# Lossy exact values remain recoverable through the local companion tag.
_HTTK_CIF_DECIMAL_DIGITS = 16


def _cif_validate_name(name_unfiltered, context=None):
    if context is not None:
        context = context + ": " + name_unfiltered
    name = _cif_validate_non_blank_char(name_unfiltered, context)
    if len(name) > 75:
        raise ValueError(f"CIF data name exceeds 75 characters: {name_unfiltered}")
    return name


def _cif_is_float(data_value):
    return _cif_float_regex.match(data_value) is not None


def _cif_is_simplestring(data_value):
    return _cif_simplestring_regex.match(data_value) is not None


def _cif_is_int(data_value):
    return _cif_integer_regex.match(data_value) is not None


def _cif_validate_non_blank_char(s, context=None):
    out = s.translate(_cif_unicode_translation_table)
    if out != s:
        if context is not None:
            warnings.warn(f"write_cif removed non-permitted characters in {context}", RuntimeWarning, stacklevel=2)
        else:
            warnings.warn("write_cif removed non-permitted characters", RuntimeWarning, stacklevel=2)
    return out


def _cif_write_semicolontextfield(f, lines, noteol, max_line_length):
    if noteol:
        f.write("\n")
        noteol = False
    for i in range(len(lines)):
        lines[i] = lines[i].rstrip("\r\n")
        if lines[i] and lines[i][0] == ';':
            warnings.warn(
                "write_cif inserted a space before a semicolon-leading text line", RuntimeWarning, stacklevel=2
            )
            lines[i] = ' ' + lines[i]
        if len(lines[i]) > max_line_length:
            f.write(";\\" + "\n")
            break
    else:
        f.write(";")
    for line in lines:
        if len(line) > max_line_length:
            sublines = [line[i : i + max_line_length - 2] for i in range(0, len(line), max_line_length - 2)]
            # Handle a wonderful corner case: the line splitting for length creates lines that start with one, or more, semi-colons..., sigh...
            for i in range(1, len(sublines)):
                if sublines[i][0] == ";":
                    if len(sublines[i]) > 1 and sublines[i][1] != ";":
                        # If its just a single semi-colon, move it to the previous line, which we saved space for by splitting at max_line_length-2
                        sublines[i - 1] += ";"
                        sublines[i] = sublines[i][1:]
                    else:
                        # Multiple semi-colons in a row, or a semi-colon + newline, this is a possibly unresolvable case (think long string of only semi-colons)
                        # fudge a solution by inserting a space
                        warnings.warn(
                            "write_cif inserted a space before a semicolon in a long text line",
                            RuntimeWarning,
                            stacklevel=2,
                        )
                        sublines[i] = " " + sublines[i]
            for subline in sublines:
                f.write(subline + "\\" + "\n")
        else:
            f.write(line + "\n")

    f.write(";\n")
    return False


def _cif_write_data_value(f, orig_data_value, noteol, max_line_length, inloop):
    if orig_data_value is None:
        f.write("?")
        return True
    else:
        data_value = str(orig_data_value)
    has_whitespace = len(data_value.split()) > 1
    lines = data_value.splitlines()
    has_lines = len(lines) > 1
    has_single_quote = data_value.find("'") != -1
    has_double_quote = data_value.find('"') != -1
    too_long = len(data_value) + 2 > max_line_length
    if has_lines or (has_single_quote and has_double_quote) or too_long:
        noteol = _cif_write_semicolontextfield(f, lines, noteol, max_line_length)
        return noteol
    elif has_double_quote or (has_whitespace and not has_single_quote) or data_value == "":
        f.write("'" + data_value + "'")
        return True
    elif has_single_quote or (has_whitespace and not has_double_quote):
        f.write('"' + data_value + '"')
        return True
    elif _cif_is_float(data_value) or _cif_is_int(data_value) or inloop and _cif_is_simplestring(data_value):
        f.write(data_value)
        return True
    else:
        f.write("'" + data_value + "'")
        return True


def is_sequence(val):
    return isinstance(val, Iterable) and not isinstance(val, str)


def _cif_write_data_block(f, data_block, max_line_length):
    for key in data_block:
        val = data_block[key]
        if key.startswith("loop_"):
            f.write("loop_\n")
            outdata_columns = []
            for unfiltered_column in val:
                column = _cif_validate_name(unfiltered_column, "column name: " + unfiltered_column)
                f.write("_" + column + "\n")
                outdata_columns += [data_block[unfiltered_column]]
            if len(outdata_columns) > 0:
                noteol = False
                for i in range(len(outdata_columns[0])):
                    column_count = 0
                    for j in range(len(outdata_columns)):
                        column_count += len(str(outdata_columns[j][i])) + 2
                        if column_count > max_line_length and noteol:
                            f.write("\n")
                            column_count = 0
                            noteol = False
                        noteol = _cif_write_data_value(f, outdata_columns[j][i], noteol, max_line_length, inloop=True)
                        if noteol:
                            f.write(" ")
                            column_count += 1
                        else:
                            column_count = 0
                    if noteol:
                        noteol = False
                        f.write("\n")
        elif is_sequence(val):
            continue
        else:
            data_name = _cif_validate_name(key)
            # Do we have space _ + key + space + quote + the whole data value + quote?, if not, preemptively break line
            f.write("_" + data_name + " ")
            if len(data_name) + len(str(val)) + 4 > max_line_length:
                f.write("\n")
                noteol = False
            else:
                noteol = True
            noteol = _cif_write_data_value(f, val, noteol, max_line_length, inloop=False)
            if noteol:
                f.write("\n")
                noteol = False


def write_cif(
    destination: str | os.PathLike[str] | io.TextIOBase,
    data,
    header: str | None = None,
    max_line_length: int = 80,
) -> None:
    """Write CIF ``data`` to a path or open text stream.

    ``data`` is an iterable of ``(block_name, block)`` pairs. A block maps data names to
    scalar values and uses ``loop_N`` keys to list its loop columns. ``header``, when
    supplied, is written before the data blocks.
    """

    with ExitStack() as stack:
        f: io.TextIOBase
        if isinstance(destination, (str, os.PathLike)):
            f = stack.enter_context(open(destination, "w", encoding="utf-8"))
        else:
            f = destination
        if header is not None:
            lines = header.splitlines()
            for line in lines:
                if len(line) > max_line_length:
                    header = "#\n" + header
                    break
            for line in lines:
                if len(line) > max_line_length:
                    sublines = [line[i : i + 79] for i in range(0, len(line), 79)]
                    for subline in sublines:
                        f.write(subline + "\\" + "\n")
                else:
                    f.write(line + "\n")

        data_block_count = -1
        for data_block in data:
            data_block_count += 1
            data_block_name_unfiltered = data_block[0]
            if data_block_name_unfiltered is None:
                data_block_name = "data_" + str(data_block_count)
            else:
                data_block_name = _cif_validate_name(data_block_name_unfiltered, "data block name")
                if data_block_name == "":
                    data_block_name = "data_" + str(data_block_count)

            f.write("data_" + data_block_name + "\n")
            _cif_write_data_block(f, data_block[1], max_line_length)


def _finite_decimal(value: Fraction) -> str | None:
    denominator = value.denominator
    while denominator % 2 == 0:
        denominator //= 2
    while denominator % 5 == 0:
        denominator //= 5
    if denominator != 1:
        return None
    places = max(
        _decimal_places(value.denominator, 2),
        _decimal_places(value.denominator, 5),
    )
    integer = value.numerator * 2 ** max(0, places - _decimal_places(value.denominator, 2))
    integer *= 5 ** max(0, places - _decimal_places(value.denominator, 5))
    sign = "-" if integer < 0 else ""
    digits = str(abs(integer)).rjust(places + 1, "0")
    if places == 0:
        return sign + digits
    return f"{sign}{digits[:-places]}.{digits[-places:]}".rstrip("0").rstrip(".")


def _decimal_places(denominator: int, factor: int) -> int:
    places = 0
    while denominator % factor == 0:
        denominator //= factor
        places += 1
    return places


def _dual_cif_value(value: object) -> tuple[object, str | None]:
    """Return a standard decimal and companion, using 16 significant digits when lossy."""
    if value is None:
        return None, None
    token = str(value).strip().strip("'\"")
    exact = Fraction(token)
    finite = _finite_decimal(exact)
    if finite is not None:
        return (finite if "/" in token or "e" in token.lower() else token), None
    with localcontext() as context:
        context.prec = _HTTK_CIF_DECIMAL_DIGITS
        decimal = Decimal(exact.numerator) / Decimal(exact.denominator)
        standard = format(decimal, f".{_HTTK_CIF_DECIMAL_DIGITS}g")
    return standard, token


def _dual_cif_column(
    values: Iterable[object], *, pad_coordinates: bool = False
) -> tuple[list[object], list[str | None], bool]:
    standard: list[object] = []
    exact: list[str | None] = []
    for value in values:
        visible, companion = _dual_cif_value(value)
        if pad_coordinates and visible is not None and companion is None:
            visible = _pad_coordinate_decimal(visible)
        standard.append(visible)
        exact.append(companion)
    return standard, exact, any(value is not None for value in exact)


def _pad_coordinate_decimal(value: object) -> object:
    """Give standard coordinate decimals a fixed 16-place precision claim."""
    text = str(value)
    whole, dot, fraction = text.partition(".")
    if not dot:
        return text + "." + "0" * _HTTK_CIF_DECIMAL_DIGITS
    return whole + "." + fraction.ljust(_HTTK_CIF_DECIMAL_DIGITS, "0")


def _neutral_cif_block(block: Mapping[str, object]) -> dict[str, object]:
    """Turn one ``read_cif_asus`` block into the low-level writer's block shape."""
    raw: dict[str, object] = {}
    cell_tags = (
        "cell_length_a",
        "cell_length_b",
        "cell_length_c",
        "cell_angle_alpha",
        "cell_angle_beta",
        "cell_angle_gamma",
    )
    exact_cell = cast(Iterable[object] | None, block.get("cell_parameters_exact"))
    if exact_cell is not None:
        for tag, value in zip(cell_tags, exact_cell):
            raw[tag], exact = _dual_cif_value(value)
            if exact is not None:
                raw[f"httk_{tag}_exact"] = exact
    for source, target in (
        ("space_group_nbr", "space_group_IT_number"),
        ("space_group_name_hm", "space_group_name_H-M_alt"),
        ("space_group_name_hall", "space_group_name_Hall"),
    ):
        if source in block:
            raw[target] = block[source]

    symops = cast(Iterable[str], block["symops_xyz"])
    raw["loop_symops"] = ["space_group_symop_operation_xyz"]
    raw["space_group_symop_operation_xyz"] = list(symops)

    positions = list(cast(Iterable[Sequence[object]], block["positions_exact"]))
    symbols = list(cast(Iterable[str], block.get("symbols", ())))
    labels = list(cast(Iterable[str], block.get("labels", symbols)))
    loop_atoms = [
        "atom_site_label",
        "atom_site_type_symbol",
        "atom_site_fract_x",
        "atom_site_fract_y",
        "atom_site_fract_z",
    ]
    raw["loop_atoms"] = loop_atoms
    raw["atom_site_label"] = labels
    raw["atom_site_type_symbol"] = symbols
    for index, tag in enumerate(("atom_site_fract_x", "atom_site_fract_y", "atom_site_fract_z")):
        values, exact_values, has_exact = _dual_cif_column((row[index] for row in positions), pad_coordinates=True)
        raw[tag] = values
        if has_exact:
            raw[f"httk_{tag}_exact"] = exact_values
            loop_atoms.append(f"httk_{tag}_exact")
    occupancies_exact = cast(Iterable[object] | None, block.get("occupancies_exact"))
    occupancies = cast(Iterable[object] | None, block.get("occupancies"))
    occupancy_values = occupancies_exact if occupancies_exact is not None else occupancies
    if occupancy_values is not None:
        values, exact_values, has_exact = _dual_cif_column(occupancy_values)
        loop_atoms.append("atom_site_occupancy")
        raw["atom_site_occupancy"] = values
        if has_exact:
            loop_atoms.append("httk_atom_site_occupancy_exact")
            raw["httk_atom_site_occupancy_exact"] = exact_values
    return raw


def _write_cif_payload(destination, data: Mapping[str, object], **kwargs: object) -> None:
    """Write the neutral CIF payload returned by ``read_cif_asus``."""
    blocks = data.get("blocks")
    if blocks is None:
        blocks = [data]
    options: dict[str, Any] = {"header": cast(str | None, data.get("header"))}
    options.update(kwargs)
    write_cif(
        destination,
        [("structure", _neutral_cif_block(block)) for block in cast(Iterable[Mapping[str, object]], blocks)],
        **options,
    )
