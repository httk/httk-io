"""Minimal example: read a CIF file through the httk-core loader dispatch.

Importing :mod:`httk.core` triggers handler discovery, which imports
``httk.handlers.io`` and registers the ``.cif`` loader.  ``httk.core.load``
then dispatches on the file extension and returns the raw parsed CIF data
(a ``(data_blocks, header)`` tuple), not a Structure object.
"""

import tempfile
from pathlib import Path

import httk.core
import httk.core.register

CIF_TEXT = """#Example NaCl-like cell
data_example
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na 0.0 0.0 0.0
Cl 0.5 0.5 0.5
"""


def main() -> None:
    print("Registered loader extensions:", httk.core.register.known_extensions())

    with tempfile.TemporaryDirectory() as tmp:
        cif_path = Path(tmp) / "example.cif"
        cif_path.write_text(CIF_TEXT, encoding="utf-8")

        data_blocks, header = httk.core.load(str(cif_path))

    print("Number of data blocks:", len(data_blocks))
    name, block = data_blocks[0]
    print("First data block name:", name)
    print("Parsed data keys:", list(block.keys()))
    print("Cell length a:", block.get("cell_length_a"))
    print("Atom site labels:", block.get("atom_site_label"))


if __name__ == "__main__":
    main()
