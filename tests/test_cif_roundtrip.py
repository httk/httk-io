"""Unit tests for the low-level cif reader/writer and numeric parsers."""

import pytest

from httk.io.cif.cif_parser import parse_cif_float, parse_cif_fraction, parse_cif_int
from httk.io.cif.cif_reader import read_cif
from httk.io.cif.cif_writer import write_cif

CIF_TEXT = """#round-trip header
data_example
_cell_length_a 4.0
_symmetry_space_group_name_h-m 'F m -3 m'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na 0.0 0.0 0.0
Cl 0.5 0.5 0.5
"""


def test_read_cif_structure(tmp_path):
    src = tmp_path / "in.cif"
    src.write_text(CIF_TEXT, encoding="utf-8")

    data, header = read_cif(str(src))

    assert header.startswith("#round-trip header")
    name, block = data[0]
    assert name == "example"
    assert block["loop_0"] == [
        "atom_site_label",
        "atom_site_fract_x",
        "atom_site_fract_y",
        "atom_site_fract_z",
    ]
    assert block["atom_site_label"] == ["Na", "Cl"]


def test_read_write_read_roundtrip(tmp_path):
    src = tmp_path / "in.cif"
    src.write_text(CIF_TEXT, encoding="utf-8")

    data, header = read_cif(str(src))

    out = tmp_path / "out.cif"
    write_cif(str(out), data, header)

    data2, _ = read_cif(str(out))

    assert data2[0][0] == data[0][0]
    assert dict(data2[0][1]) == dict(data[0][1])


def test_writer_uses_question_mark_for_none_and_reader_restores_it(tmp_path):
    out = tmp_path / "none.cif"
    write_cif(str(out), [("example", {"cell_length_a": None})])

    _name, block = read_cif(str(out))[0][0]
    assert block["cell_length_a"] == "?"
    assert parse_cif_float(block["cell_length_a"]) is None


def test_writer_handles_blank_lines_in_multiline_values(tmp_path):
    out = tmp_path / "blank-line.cif"
    write_cif(str(out), [("example", {"audit_creation_method": "first\n\nlast"})])

    _name, block = read_cif(str(out))[0][0]
    assert block["audit_creation_method"] == "first\n\nlast"


def test_writer_rejects_overlong_data_names(tmp_path):
    with pytest.raises(ValueError, match="75"):
        write_cif(str(tmp_path / "long.cif"), [("example", {"a" * 76: "value"})])


def test_parse_cif_float_with_uncertainty():
    assert parse_cif_float("1.234(5)") == 1.234
    assert parse_cif_float("-0.5") == -0.5
    assert parse_cif_float("?") is None
    assert abs(parse_cif_float("1/3") - (1.0 / 3.0)) < 1e-9


def test_parse_cif_fraction_preserves_fraction_and_decimal_tokens():
    from fractions import Fraction

    assert parse_cif_fraction("1/3") == Fraction(1, 3)
    assert parse_cif_fraction("0.125") == Fraction(1, 8)


def test_malformed_numeric_salvage_requires_pragmatic_mode():
    with pytest.raises(ValueError, match="Invalid CIF numeric"):
        parse_cif_float("1.2(bad)")
    with pytest.warns(RuntimeWarning, match="Salvaging"):
        assert parse_cif_float("1.2(bad)", pragmatic=True) == 1.2


def test_parse_cif_float_meta_reports_esd():
    value, meta = parse_cif_float("1.234(5)", meta=True)
    assert value == 1.234
    assert meta["esd"] is not None
    assert meta["esd"] > 0.0


def test_parse_cif_int():
    assert parse_cif_int("123(4)") == 123
    assert parse_cif_int("3E2") == 300
