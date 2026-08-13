import logging
from pathlib import Path

import pytest
from httk.core._plugins import resolve_callable

from httk.io.cif.cif_parser import read_cif_asus
from httk.io.cif.cif_reader import _PROTECTED_LOOP_TAGS, read_cif

FIXTURE = Path(__file__).parent / "fixtures" / "malformed_auxiliary_loop.cif"
HINT = " (an auxiliary loop like this can be dropped by loading with autocorrect=True, which applies documented repairs with warnings)"


def test_protected_set_includes_atom_declaration_tags():
    assert {
        "atom_site_wyckoff_label",
        "atom_site_symmetry_multiplicity",
    } <= _PROTECTED_LOOP_TAGS


def test_autocorrect_drops_malformed_auxiliary_loop_and_stamps_payload(caplog):
    with pytest.raises(ValueError) as error:
        read_cif(FIXTURE)
    assert HINT in str(error.value)

    with caplog.at_level(logging.WARNING):
        reader = resolve_callable("httk.io.cif:read_cif_asus")
        assert reader is read_cif_asus
        payload = reader(FIXTURE, autocorrect=True)

    assert payload["autocorrect"] is True
    assert payload["blocks"][0]["labels"] == ["Na1"]
    assert len(caplog.records) == 1
    record = caplog.records[0]
    assert record.name == "httk.io.cif.cif_reader"
    assert record.context == "cif"
    assert "_audit_tag" in record.message
    assert "dropped" in record.message


def test_autocorrect_does_not_drop_malformed_atom_site_loop(tmp_path):
    atom_site_fixture = tmp_path / "malformed_atom_site.cif"
    atom_site_fixture.write_text(FIXTURE.read_text(encoding="utf-8").replace("_audit_tag", "_atom_site_occupancy"))

    for autocorrect in (False, True):
        with pytest.raises(ValueError) as error:
            read_cif(atom_site_fixture, autocorrect=autocorrect)
        assert HINT not in str(error.value)


def test_autocorrect_does_not_drop_malformed_modulation_loop(tmp_path):
    modulation_fixture = tmp_path / "malformed_modulation.cif"
    modulation_fixture.write_text(FIXTURE.read_text(encoding="utf-8").replace("_audit_tag", "_cell_wave_vector_x"))

    for autocorrect in (False, True):
        with pytest.raises(ValueError) as error:
            read_cif(modulation_fixture, autocorrect=autocorrect)
        assert HINT not in str(error.value)


def test_autocorrect_does_not_drop_malformed_wyckoff_declaration_loop(tmp_path):
    wyckoff_fixture = tmp_path / "malformed_wyckoff.cif"
    wyckoff_fixture.write_text(FIXTURE.read_text(encoding="utf-8").replace("_audit_tag", "_atom_site_Wyckoff_label"))

    for autocorrect in (False, True):
        with pytest.raises(ValueError) as error:
            read_cif(wyckoff_fixture, autocorrect=autocorrect)
        assert HINT not in str(error.value)
