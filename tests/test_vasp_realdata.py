"""Workspace-only checks against the external VASP example-data tree."""

from pathlib import Path

import pytest

from httk.io.vasp import OutcarFile, VASPOutputs, XdatcarFile, read_oszicar, read_poscar

FIXTURES = Path(__file__).resolve().parent.parent.parent / "electronic-structure-example-data"
pytestmark = pytest.mark.skipif(not FIXTURES.exists(), reason="workspace-only real-data fixtures not present")


def test_static_nacl_outcar() -> None:
    outcar = OutcarFile(FIXTURES / "Static/VASP/NaCl/OUTCAR.bz2")
    assert outcar.version_string.startswith("vasp.")
    assert outcar.xc == "Perdew-Burke-Ernzerhof (PBE)"
    assert outcar.final_energies.free_energy is not None
    assert outcar.final_energies.energy_without_entropy is not None
    assert outcar.final_energies.energy_sigma0 is not None
    assert outcar.final_energies.final is True
    assert outcar.completed is True


def test_partial_charges_build_suffix_version() -> None:
    outcar = OutcarFile(FIXTURES / "PartialCharges/VASP/C-diamond/OUTCAR.bz2")
    assert outcar.version_string.startswith("vasp.5.4.4.18Apr17")


def test_static_cu_fcc_version_numbers() -> None:
    outcar = OutcarFile(FIXTURES / "Static/VASP/Cu-FCC/OUTCAR.bz2")
    assert outcar.version_numbers[:3] == (5, 2, 12)


@pytest.mark.extended
def test_md_outcar_has_10000_frames() -> None:
    outcar = OutcarFile(FIXTURES / "MD/VASP/Al_300K/OUTCAR.bz2")
    assert outcar.nframes == 10000


def test_md_xdatcar_oszicar_and_poscar_agree() -> None:
    directory = FIXTURES / "MD/VASP/Al_300K"
    xdatcar = XdatcarFile(directory / "XDATCAR.bz2")
    assert xdatcar.nframes == 10000
    oszicar = read_oszicar(directory / "OSZICAR.bz2")
    assert len(oszicar["ionic_steps"]) == xdatcar.nframes
    poscar = read_poscar(directory / "POSCAR")
    first = next(xdatcar.frames())
    assert len(first["coords"]) == sum(poscar["counts"])


def test_md_vasp_outputs_composite() -> None:
    outputs = VASPOutputs(FIXTURES / "MD/VASP/Al_300K")
    assert outputs.poscar is not None
    assert outputs.outcar is not None
    assert outputs.xdatcar is not None
    assert outputs.oszicar is not None
    assert outputs.potcar is None
    outputs.close()
