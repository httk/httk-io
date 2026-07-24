"""Tests for the VASP POSCAR/CONTCAR reader and its loader registration."""

import bz2
from pathlib import Path

import pytest

import httk.core
from httk.io.vasp import read_poscar

VASP5_SELECTIVE = """SmFeO3 slab
1.0
1.0 0.0 0.0
0.0 2.0 0.0
0.0 0.0 3.0
Si O
1 2
Selective dynamics
Direct
0.0 0.0 0.0 T T F
0.25 0.25 0.25 F F F
0.5 0.5 0.5 T T T
0.0000000000000000  0.0000000000000000  0.0000000000000000
0.1 0.1 0.1
"""

VASP4_CARTESIAN_LABELS = """VASP4 cartesian with trailing labels
2.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
2 1
Cartesian
0.0 0.0 0.0 Si
0.5 0.5 0.5 Si
0.25 0.25 0.25 O
"""

NEGATIVE_SCALE = """volume-scaled cell
-16.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
He
1
Direct
0.0 0.0 0.0
"""


def test_vasp5_selective_dynamics_and_velocity_block() -> None:
    data = read_poscar(VASP5_SELECTIVE.splitlines(keepends=True))
    assert data["format"] == "vasp-poscar"
    assert data["comment"] == "SmFeO3 slab"
    assert data["scale"] == "1.0"
    assert data["volume"] is None
    assert data["symbols"] == ["Si", "O"]
    assert data["counts"] == [1, 2]
    assert data["cartesian"] is False
    assert data["coords"] == [["0.0", "0.0", "0.0"], ["0.25", "0.25", "0.25"], ["0.5", "0.5", "0.5"]]
    # The trailing velocity block is ignored (only 3 coordinate rows read).
    assert data["selective_dynamics"] == [[True, True, False], [False, False, False], [True, True, True]]


def test_vasp4_cartesian_trailing_labels() -> None:
    data = read_poscar(VASP4_CARTESIAN_LABELS.splitlines(keepends=True))
    assert data["symbols"] is None  # VASP-4: no species line
    assert data["counts"] == [2, 1]
    assert data["cartesian"] is True
    assert data["selective_dynamics"] is None
    # Trailing per-line species labels are dropped (first three tokens only).
    assert data["coords"] == [["0.0", "0.0", "0.0"], ["0.5", "0.5", "0.5"], ["0.25", "0.25", "0.25"]]


def test_negative_scale_is_volume() -> None:
    data = read_poscar(NEGATIVE_SCALE.splitlines(keepends=True))
    assert data["scale"] is None
    assert data["volume"] == "16.0"


def test_malformed_reports_line_number() -> None:
    broken = "c\n1.0\n1 0 0\n0 1 0\n"  # truncated before the third cell row
    with pytest.raises(ValueError) as excinfo:
        read_poscar(broken.splitlines(keepends=True))
    assert "line 5" in str(excinfo.value)


# --- loader registration + end-to-end load() ----------------------------------


def test_poscar_registration() -> None:
    assert ".poscar" in httk.core.register.known_extensions()
    assert ".vasp" in httk.core.register.known_extensions()
    assert "poscar" in httk.core.register.known_filenames()
    assert "contcar" in httk.core.register.known_filenames()


def test_load_contcar_by_basename(tmp_path: Path) -> None:
    contcar = tmp_path / "CONTCAR"
    contcar.write_text(VASP5_SELECTIVE, encoding="utf-8")
    data = httk.core.load(str(contcar))
    assert data["format"] == "vasp-poscar"
    assert data["counts"] == [1, 2]


def test_load_contcar_bz2_transparent(tmp_path: Path) -> None:
    contcar_bz2 = tmp_path / "CONTCAR.bz2"
    contcar_bz2.write_bytes(bz2.compress(VASP5_SELECTIVE.encode("utf-8")))
    data = httk.core.load(str(contcar_bz2))
    assert data["format"] == "vasp-poscar"
    assert data["symbols"] == ["Si", "O"]


CIF_TEXT = """#header
data_x
_cell_length_a 5.64
"""


def test_load_cif_bz2_transparent(tmp_path: Path) -> None:
    cif_bz2 = tmp_path / "sample.cif.bz2"
    cif_bz2.write_bytes(bz2.compress(CIF_TEXT.encode("utf-8")))
    datalist, header = httk.core.load(str(cif_bz2))
    assert header.startswith("#header")
    name, block = datalist[0]
    assert block["cell_length_a"] == "5.64"
