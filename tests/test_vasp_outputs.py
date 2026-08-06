"""Tests for the lazy VASP directory composite."""

import bz2
from pathlib import Path

import httk.core
import pytest

from httk.io.vasp import OutcarFile, VASPOutputs

POSCAR = """synthetic POSCAR
1.0
1 0 0
0 1 0
0 0 1
Si
1
Direct
0 0 0
"""
OSZICAR = "  1 F= -.10000 E0= -.20000 d E =-.00100\n"
OUTCAR = " vasp.5.2.12 synthetic\n General timing and accounting informations\n"


def test_outputs_lazy_resolution_and_suffixes(tmp_path: Path) -> None:
    (tmp_path / "POSCAR").write_text(POSCAR, encoding="utf-8")
    (tmp_path / "CONTCAR.bz2").write_bytes(bz2.compress(POSCAR.encode("utf-8")))
    (tmp_path / "OUTCAR.bz2").write_bytes(bz2.compress(OUTCAR.encode("utf-8")))
    (tmp_path / "OSZICAR").write_text(OSZICAR, encoding="utf-8")
    outputs = VASPOutputs(tmp_path)
    assert outputs.poscar is not None
    assert outputs.contcar is not None
    assert isinstance(outputs.outcar, OutcarFile)
    assert outputs.xdatcar is None
    assert outputs.oszicar is not None
    assert outputs.potcar is None
    loaded = httk.core.load(str(tmp_path / "OSZICAR"), raw=True)
    assert loaded["format"] == "vasp-oszicar"
    assert not outputs.closed


def test_outputs_owns_file_objects_and_close_is_idempotent(tmp_path: Path) -> None:
    (tmp_path / "OUTCAR").write_text(OUTCAR, encoding="utf-8")
    with VASPOutputs(tmp_path) as outputs:
        outcar = outputs.outcar
        assert outcar is not None
    assert outputs.closed
    assert outcar is not None and outcar.closed
    outputs.close()
    with pytest.raises(ValueError, match="closed VASP outputs"):
        _ = outputs.outcar
