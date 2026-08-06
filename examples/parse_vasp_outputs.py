"""Read a small synthetic VASP output directory."""

from pathlib import Path
from tempfile import TemporaryDirectory

from httk.io.vasp import VASPOutputs

POSCAR = """synthetic
1.0
1 0 0
0 1 0
0 0 1
He
1
Direct
0 0 0
"""
OSZICAR = "  1 F= -.1 E0= -.2 dE =-.01\n"
XDATCAR = "\n".join(POSCAR.splitlines()[:7]) + "\nDirect configuration= 1\n0 0 0\n"
OUTCAR = """ vasp.5.2.12 synthetic
   GGA     =    PE
   ENCUT  = 300.0 eV
 ----------------------------------------- Iteration    1(   1)  ---------------------------------------
   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
   free  energy   TOTEN  =       -0.100 eV
   energy  without entropy=      -0.100  energy(sigma->0) =      -0.100
 General timing and accounting informations
"""


def main() -> None:
    with TemporaryDirectory() as temporary:
        directory = Path(temporary)
        (directory / "POSCAR").write_text(POSCAR, encoding="utf-8")
        (directory / "OSZICAR").write_text(OSZICAR, encoding="utf-8")
        (directory / "XDATCAR").write_text(XDATCAR, encoding="utf-8")
        (directory / "OUTCAR").write_text(OUTCAR, encoding="utf-8")

        with VASPOutputs(directory) as outputs:
            assert outputs.poscar is not None
            assert outputs.oszicar is not None
            assert outputs.xdatcar is not None
            assert outputs.outcar is not None
            assert outputs.xdatcar.nframes == 1
            assert outputs.outcar.final_energies.final
            print("VASP outputs:", outputs.poscar["symbols"], outputs.xdatcar.nframes)


if __name__ == "__main__":
    main()
