"""Unit tests for the low-level cif reader/writer and numeric parsers."""

from httk.io.cif.cif_parser import parse_cif_float, parse_cif_int
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


def test_parse_cif_float_with_uncertainty():
    assert parse_cif_float("1.234(5)") == 1.234
    assert parse_cif_float("-0.5") == -0.5
    assert parse_cif_float("?") is None
    assert abs(parse_cif_float("1/3") - (1.0 / 3.0)) < 1e-9


def test_parse_cif_float_meta_reports_esd():
    value, meta = parse_cif_float("1.234(5)", meta=True)
    assert value == 1.234
    assert meta["esd"] is not None
    assert meta["esd"] > 0.0


def test_parse_cif_int():
    assert parse_cif_int("123(4)") == 123
    assert parse_cif_int("3E2") == 300
