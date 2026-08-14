"""End-to-end: httk.core.load dispatches .cif to the httk-io reader."""

from pathlib import Path

import httk.core

CIF_TEXT = """#a small header
data_nacl
_cell_length_a 5.64
_cell_length_b 5.64
_cell_length_c 5.64
_symmetry_space_group_name_h-m 'F m -3 m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na1 Na 0.0 0.0 0.0
Cl1 Cl 0.5 0.5 0.5
"""


def _write_cif(tmp_path: Path) -> Path:
    cif_path = tmp_path / "nacl.cif"
    cif_path.write_text(CIF_TEXT, encoding="utf-8")
    return cif_path


def test_load_returns_a_tagged_cif_payload(tmp_path):
    payload = httk.core.load(str(_write_cif(tmp_path)), raw=True)

    assert payload["format"] == "cif"
    assert payload["header"].startswith("#a small header")


def test_load_tolerates_a_block_that_is_not_a_structure(tmp_path):
    """CIF is a general-purpose format, so loading must not insist on crystallography.

    This block names atom sites but gives no cell and no symmetry operations, so it is not
    a structure. Loading still succeeds and records why it could not be interpreted,
    rather than failing the whole file.
    """
    payload = httk.core.load(str(_write_cif(tmp_path)), raw=True)

    assert payload["blocks"] == []
    assert [item["block"] for item in payload["unparsed"]] == ["nacl"]
    assert payload["unparsed"][0]["reason"]


def test_the_low_level_tokenizer_is_still_available(tmp_path):
    from httk.io.cif import read_cif

    data_blocks, header = read_cif(str(_write_cif(tmp_path)))

    assert header.startswith("#a small header")
    assert len(data_blocks) == 1
    name, block = data_blocks[0]
    assert name == "nacl"
    assert block["cell_length_a"] == "5.64"
    assert block["atom_site_label"] == ["Na1", "Cl1"]
    assert block["atom_site_type_symbol"] == ["Na", "Cl"]
    assert block["atom_site_fract_x"] == ["0.0", "0.5"]


def test_cif_payload_preserves_raw_symmetry_operation_tokens(tmp_path):
    path = tmp_path / "raw.cif"
    path.write_text(
        "data_raw\n"
        "_cell_length_a 1\n_cell_length_b 1\n_cell_length_c 1\n"
        "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
        "loop_\n_space_group_symop_operation_xyz\n' x , y , z '\n"
        "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n"
        "_atom_site_fract_y\n_atom_site_fract_z\nSi1 Si 0 0 0\n",
        encoding="utf-8",
    )
    block = httk.core.load(path, raw=True)["blocks"][0]
    assert block["symops_xyz"] == (" x , y , z ",)
    assert "symops" not in block


def test_load_unknown_extension_raises(tmp_path):
    bad = tmp_path / "data.unknownext"
    bad.write_text("nonsense", encoding="utf-8")
    try:
        httk.core.load(str(bad))
    except Exception:
        return
    raise AssertionError("expected load() to reject an unregistered extension")


def test_load_mcif_returns_a_tagged_neutral_payload(tmp_path):
    mcif = tmp_path / "magnetic.mcif"
    mcif.write_text(
        (Path(__file__).with_name("fixtures") / "magnetic_centered.mcif").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    # raw=True keeps the assertion valid when a domain module (httk-atomistic)
    # is co-installed and its format adapter would otherwise convert the payload.
    payload = httk.core.load(str(mcif), raw=True)
    assert payload["format"] == "mcif"
    assert len(payload["blocks"]) == 1
    assert payload["unparsed"] == []
    block = payload["blocks"][0]
    assert block["symops_xyz"] == ("x,y,z,+1",)
    assert "basis" not in block
    assert "positions" not in block
    assert "symops" not in block


def test_missing_atom_site_columns_are_named(tmp_path):
    path = tmp_path / "missing-symbol.cif"
    path.write_text(
        "data_missing\n"
        "_cell_length_a 1\n_cell_length_b 1\n_cell_length_c 1\n"
        "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
        "loop_\n_space_group_symop_operation_xyz\n'x,y,z'\n"
        "loop_\n_atom_site_label\n_atom_site_fract_x\n"
        "_atom_site_fract_y\n_atom_site_fract_z\nSi1 0 0 0\n",
        encoding="utf-8",
    )

    payload = httk.core.load(path, raw=True)

    assert payload["unparsed"][0]["reason"] == (
        "ValueError: CIF block is missing required atom-site column: _atom_site_type_symbol"
    )
