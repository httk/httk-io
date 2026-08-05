"""Tests for the VASP POSCAR/CONTCAR reader and its loader registration."""

import bz2
import io
from pathlib import Path

import httk.core
import pytest

from httk.io.vasp import read_poscar
from httk.io.vasp.poscar_writer import _write_poscar_payload

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


def test_poscar_writer_registration() -> None:
    for name in ("POSCAR", "CONTCAR", "x.vasp", "y.poscar", "POSCAR.bz2"):
        assert httk.core.has_writer_for(name)


def _without_precision(data: dict) -> dict:
    return {key: value for key, value in data.items() if not key.endswith("_precision")}


@pytest.mark.parametrize("source", (VASP5_SELECTIVE, VASP4_CARTESIAN_LABELS, NEGATIVE_SCALE))
def test_poscar_writer_round_trips_payload(source: str, tmp_path: Path) -> None:
    payload = read_poscar(source.splitlines(keepends=True))
    path = tmp_path / "POSCAR"
    httk.core.save(payload, path)
    assert _without_precision(read_poscar(path)) == _without_precision(payload)


def test_poscar_writer_round_trips_zero_count_species(tmp_path: Path) -> None:
    payload = read_poscar(VASP5_SELECTIVE.splitlines(keepends=True))
    payload["counts"] = [2, 0]
    payload["coords"] = payload["coords"][:2]
    payload["selective_dynamics"] = payload["selective_dynamics"][:2]
    path = tmp_path / "POSCAR"
    httk.core.save(payload, path)
    assert _without_precision(read_poscar(path)) == _without_precision(payload)


def test_poscar_writer_round_trips_compressed_destination(tmp_path: Path) -> None:
    payload = read_poscar(VASP5_SELECTIVE.splitlines(keepends=True))
    path = tmp_path / "POSCAR.bz2"
    httk.core.save(payload, path)
    assert _without_precision(httk.core.load(str(path), raw=True)) == _without_precision(payload)


def test_poscar_writer_round_trips_open_stream() -> None:
    payload = read_poscar(VASP5_SELECTIVE.splitlines(keepends=True))
    stream = io.StringIO()
    _write_poscar_payload(stream, payload)
    stream.seek(0)
    assert _without_precision(read_poscar(stream)) == _without_precision(payload)


@pytest.mark.parametrize(
    "update",
    (
        {"volume": "16.0"},
        {"scale": None, "volume": None},
        {"cell": [["1.0", "0.0", "0.0"], ["0.0", "1.0", "0.0"]]},
        {"counts": [1, 1]},
        {"symbols": ["Si"]},
        {"selective_dynamics": [[True, True, False], [False, False, False]]},
    ),
)
def test_poscar_writer_rejects_malformed_payload(update: dict, tmp_path: Path) -> None:
    payload = read_poscar(VASP5_SELECTIVE.splitlines(keepends=True))
    payload.update(update)
    with pytest.raises(ValueError):
        _write_poscar_payload(tmp_path / "POSCAR", payload)


def test_load_contcar_by_basename(tmp_path: Path) -> None:
    contcar = tmp_path / "CONTCAR"
    contcar.write_text(VASP5_SELECTIVE, encoding="utf-8")
    data = httk.core.load(str(contcar), raw=True)
    assert data["format"] == "vasp-poscar"
    assert data["counts"] == [1, 2]


def test_load_contcar_bz2_transparent(tmp_path: Path) -> None:
    contcar_bz2 = tmp_path / "CONTCAR.bz2"
    contcar_bz2.write_bytes(bz2.compress(VASP5_SELECTIVE.encode("utf-8")))
    data = httk.core.load(str(contcar_bz2), raw=True)
    assert data["format"] == "vasp-poscar"
    assert data["symbols"] == ["Si", "O"]


CIF_TEXT = """#header
data_x
_cell_length_a 5.64
"""


def test_load_cif_bz2_transparent(tmp_path: Path) -> None:
    cif_bz2 = tmp_path / "sample.cif.bz2"
    cif_bz2.write_bytes(bz2.compress(CIF_TEXT.encode("utf-8")))
    payload = httk.core.load(str(cif_bz2), raw=True)
    assert payload["format"] == "cif"
    assert payload["header"].startswith("#header")

    # The point of this test is transparent decompression, so check it at the tokenizer
    # too, where the tags are visible verbatim.
    from httk.io.cif import read_cif

    datalist, _header = read_cif(str(cif_bz2))
    _name, block = datalist[0]
    assert block["cell_length_a"] == "5.64"
