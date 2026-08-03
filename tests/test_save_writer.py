from httk.core import load, save

from httk.io.cif.cif_parser import parse_cif_float
from httk.io.cif.cif_reader import read_cif


def test_registered_cif_writer_round_trips_neutral_payload(tmp_path):
    source = {
        "format": "cif",
        "blocks": [
            {
                "cell_parameters_exact": ("5.64",) * 3 + ("90",) * 3,
                "positions_exact": (("1/3", "0", "0"),),
                "symops_xyz": ("x,y,z",),
                "symbols": ("Si",),
                "labels": ("Si1",),
                "occupancies_exact": ("1",),
                "space_group_nbr": "1",
                "space_group_name_hm": "P 1",
            }
        ],
    }
    path = tmp_path / "roundtrip.cif"
    save(source, path)
    result = load(path, raw=True)["blocks"][0]
    assert result["positions_exact"] == [("1/3", "0.0000000000000000", "0.0000000000000000")]
    assert result["cell_parameters_exact"] == ("5.64",) * 3 + ("90",) * 3
    _name, raw = read_cif(path)[0][0]
    assert raw["atom_site_fract_x"] == ["0.3333333333333333"]
    assert parse_cif_float(raw["atom_site_fract_x"][0]) == 0.3333333333333333
    assert raw["httk_atom_site_fract_x_exact"] == ["1/3"]


def test_registered_cif_writer_canonicalizes_scientific_exact_tokens(tmp_path):
    source = {
        "format": "cif",
        "blocks": [
            {
                "cell_parameters_exact": ("5.64",) * 3 + ("90",) * 3,
                "positions_exact": (("1e-3", "1E2", "-1.5e-4"),),
                "symops_xyz": ("x,y,z",),
                "symbols": ("Si",),
                "labels": ("Si1",),
                "space_group_nbr": "1",
                "space_group_name_hm": "P 1",
            }
        ],
    }
    path = tmp_path / "scientific.cif"
    save(source, path)
    _name, raw = read_cif(path)[0][0]
    assert raw["atom_site_fract_x"] == ["0.0010000000000000"]
    assert raw["atom_site_fract_y"] == ["100.0000000000000000"]
    assert raw["atom_site_fract_z"] == ["-0.0001500000000000"]
    assert all(
        parse_cif_float(raw[tag][0]) is not None
        for tag in ("atom_site_fract_x", "atom_site_fract_y", "atom_site_fract_z")
    )
