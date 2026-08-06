"""Tests for the VASP POTCAR header reader."""

import bz2
from pathlib import Path

import httk.core

from httk.io.vasp import read_potcar_summary

POTCAR_SUMMARY = """   TITEL  = PAW_PBE Na_pv 05Jan2001
   LEXCH  = PE
   POMASS = 22.990; ZVAL = 7.000
   ENMAX  = 259.561; ENMIN = 194.671 eV
   licensed body omitted
   TITEL  = PAW_PBE Cl 17Jan2003
   POMASS = 35.453; ZVAL = 7.000
   ENMAX  = 280.000; ENMIN = 196.854 eV
"""


def test_potcar_summary_extracts_headers_only() -> None:
    data = read_potcar_summary(POTCAR_SUMMARY.splitlines(keepends=True))
    assert data == {
        "format": "vasp-potcar",
        "potentials": [
            {
                "titel": "PAW_PBE Na_pv 05Jan2001",
                "symbol": "Na",
                "zval": "7.000",
                "pomass": "22.990",
                "enmax": "259.561",
                "lexch": "PE",
            },
            {
                "titel": "PAW_PBE Cl 17Jan2003",
                "symbol": "Cl",
                "zval": "7.000",
                "pomass": "35.453",
                "enmax": "280.000",
                "lexch": None,
            },
        ],
    }
    assert "licensed body omitted" not in repr(data)


def test_potcar_summary_registration_and_bz2(tmp_path: Path) -> None:
    source = tmp_path / "POTCAR.summary.bz2"
    source.write_bytes(bz2.compress(POTCAR_SUMMARY.encode("utf-8")))
    data = httk.core.load(str(source), raw=True)
    assert data["format"] == "vasp-potcar"
    assert data["potentials"][0]["symbol"] == "Na"
    assert ".potcar" in httk.core.register.known_extensions()
    assert "potcar.summary" in httk.core.register.known_filenames()
