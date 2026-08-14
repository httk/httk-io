"""Inferring atom-site element symbols from labels when a CIF omits the type-symbol column."""

from pathlib import Path

import pytest

from httk.io.cif.cif_parser import _symbol_from_label, read_cif_asus, single_asu_from_cif_file

# COD diopside; its atom_site loop carries no _atom_site_type_symbol column.
COD_DIOPSIDE = Path(__file__).resolve().parents[2] / "DATA/COD/cif/1/00/00/1000008.cif"

_CELL = """data_example
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
"""


def _loop(labels):
    lines = [
        "loop_",
        "_atom_site_label",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
    ]
    for index, label in enumerate(labels):
        coord = 0.1 * index
        lines.append(f"{label} {coord} {coord} {coord}")
    return _CELL + "\n".join(lines) + "\n"


@pytest.mark.parametrize(
    ("label", "expected"),
    [
        ("MgM1", "Mg"),
        ("CaM2", "Ca"),
        ("SiT", "Si"),
        ("O1", "O"),
        ("OW1", "O"),
        ("OSi1", "O"),
        ("MGM1", "Mg"),
        ("Fe3+", "Fe"),
        ("xyz", None),
        ("1ab", None),
    ],
)
def test_symbol_from_label(label, expected):
    assert _symbol_from_label(label) == expected


def test_symbols_inferred_from_labels_with_warning(tmp_path):
    src = tmp_path / "no_type_symbol.cif"
    src.write_text(_loop(["MgM1", "CaM2", "SiT", "O1", "OW1"]), encoding="utf-8")

    with pytest.warns(RuntimeWarning, match="inferred from _atom_site_label"):
        payload = read_cif_asus(str(src))

    assert payload["blocks"], payload["unparsed"]
    assert payload["blocks"][0]["symbols"] == ["Mg", "Ca", "Si", "O", "O"]


def test_uninferable_label_raises_naming_it(tmp_path):
    src = tmp_path / "bad_label.cif"
    src.write_text(_loop(["MgM1", "Zz1"]), encoding="utf-8")

    with pytest.raises(ValueError, match="Zz1"):
        single_asu_from_cif_file(str(src))


def test_present_type_symbol_column_is_unchanged(tmp_path, recwarn):
    text = _CELL + (
        "loop_\n"
        "_atom_site_type_symbol\n"
        "_atom_site_label\n"
        "_atom_site_fract_x\n"
        "_atom_site_fract_y\n"
        "_atom_site_fract_z\n"
        "Mg Q1 0.0 0.0 0.0\n"
        "O Zz2 0.5 0.5 0.5\n"
    )
    src = tmp_path / "with_type_symbol.cif"
    src.write_text(text, encoding="utf-8")

    payload = read_cif_asus(str(src))

    assert payload["blocks"][0]["symbols"] == ["Mg", "O"]
    assert not [w for w in recwarn.list if "inferred from _atom_site_label" in str(w.message)]


@pytest.mark.skipif(not COD_DIOPSIDE.exists(), reason="workspace-only real-data fixture not present")
def test_cod_diopside_reads_without_type_symbol_column():
    with pytest.warns(RuntimeWarning, match="inferred from _atom_site_label"):
        payload = read_cif_asus(str(COD_DIOPSIDE))

    assert payload["blocks"], payload["unparsed"]
    assert payload["blocks"][0]["symbols"] == ["Mg", "Ca", "Si", "O", "O", "O"]
