"""Tests for the VASP OSZICAR reader."""

import bz2
from pathlib import Path

import httk.core

from httk.io.vasp import read_oszicar

OSZICAR = """       N       E                     dE             d eps       ncg     rms          rms(c)
DAV:   1     -1.000E+00    -1.000E+00   -2.000E+00   10   0.100E+00    0.010E+00
RMM:   2     -1.100E+00    -0.100E+00   -0.020E+00   11   0.020E+00
  1 F= -.900E+00 E0= -.800E+00 d E =-.100E+00 mag= 1.25
CG:    1     -2.000E+00     0.100E+00   -0.003E+00   12   0.030E+00
  2 F= -.950E+00 E0= -.850E+00 dE =-.050E+00
"""


def test_oszicar_groups_electronic_lines_with_following_ionic_step() -> None:
    data = read_oszicar(OSZICAR.splitlines(keepends=True))
    assert data == {
        "format": "vasp-oszicar",
        "ionic_steps": [
            {
                "n": 1,
                "F": "-.900E+00",
                "E0": "-.800E+00",
                "dE": "-.100E+00",
                "mag": "1.25",
                "electronic": [
                    {
                        "scheme": "DAV",
                        "n": 1,
                        "E": "-1.000E+00",
                        "dE": "-1.000E+00",
                        "d_eps": "-2.000E+00",
                        "ncg": "10",
                        "rms": "0.100E+00",
                        "rms_c": "0.010E+00",
                    },
                    {
                        "scheme": "RMM",
                        "n": 2,
                        "E": "-1.100E+00",
                        "dE": "-0.100E+00",
                        "d_eps": "-0.020E+00",
                        "ncg": "11",
                        "rms": "0.020E+00",
                        "rms_c": None,
                    },
                ],
            },
            {
                "n": 2,
                "F": "-.950E+00",
                "E0": "-.850E+00",
                "dE": "-.050E+00",
                "mag": None,
                "electronic": [
                    {
                        "scheme": "CG",
                        "n": 1,
                        "E": "-2.000E+00",
                        "dE": "0.100E+00",
                        "d_eps": "-0.003E+00",
                        "ncg": "12",
                        "rms": "0.030E+00",
                        "rms_c": None,
                    }
                ],
            },
        ],
        "issues": [],
    }


def test_oszicar_trailing_electronic_block_and_issues() -> None:
    text = OSZICAR[: OSZICAR.rfind("  2 F=")] + "DAV: malformed\nnot an OSZICAR line\n"
    data = read_oszicar(text.splitlines(keepends=True))
    assert data["ionic_steps"][-1] == {
        "n": None,
        "F": None,
        "E0": None,
        "dE": None,
        "mag": None,
        "electronic": [
            {
                "scheme": "CG",
                "n": 1,
                "E": "-2.000E+00",
                "dE": "0.100E+00",
                "d_eps": "-0.003E+00",
                "ncg": "12",
                "rms": "0.030E+00",
                "rms_c": None,
            },
        ],
    }
    assert any("line 6" in issue and "malformed electronic" in issue for issue in data["issues"])
    assert any("line 7" in issue and "unrecognized" in issue for issue in data["issues"])


def test_oszicar_malformed_boundary_drops_pending_iterations() -> None:
    text = """DAV:   1 -1.000E+00 -1.000E+00 -2.000E+00 10 0.100E+00
  1 E0= -.800E+00 dE =-.100E+00
DAV:   2 -2.000E+00 -2.000E+00 -3.000E+00 11 0.200E+00
  2 F= -.950E+00 E0= -.850E+00 dE =-.050E+00
"""
    data = read_oszicar(text.splitlines(keepends=True))
    assert len(data["ionic_steps"]) == 1
    assert [line["n"] for line in data["ionic_steps"][0]["electronic"]] == [2]
    assert any("dropped electronic iterations from lines 1-1" in issue for issue in data["issues"])


def test_oszicar_registration_and_bz2(tmp_path: Path) -> None:
    source = tmp_path / "OSZICAR.bz2"
    source.write_bytes(bz2.compress(OSZICAR.encode("utf-8")))
    data = httk.core.load(str(source), raw=True)
    assert data["format"] == "vasp-oszicar"
    assert data["ionic_steps"][1]["F"] == "-.950E+00"
    assert ".oszicar" in httk.core.register.known_extensions()
    assert "oszicar" in httk.core.register.known_filenames()
