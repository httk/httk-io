"""Minimal example: read a VASP POSCAR/CONTCAR through the httk-core loader.

Importing :mod:`httk.core` triggers handler discovery, which imports
``httk.handlers.io`` and registers the POSCAR loader under the ``.poscar`` /
``.vasp`` extensions and the ``POSCAR`` / ``CONTCAR`` basenames.
``httk.core.load`` then dispatches the file (transparently decompressing a
``.bz2`` file) and returns the neutral, string-preserving mapping.
"""

import bz2
import tempfile
from pathlib import Path

import httk.core
import httk.core.register

CONTCAR_TEXT = """Si primitive cell
1.0
0.0 2.7 2.7
2.7 0.0 2.7
2.7 2.7 0.0
Si
2
Direct
0.00 0.00 0.00
0.25 0.25 0.25
"""


def main() -> None:
    print("Registered loader extensions:", httk.core.register.known_extensions())
    print("Registered loader filenames:", httk.core.register.known_filenames())

    with tempfile.TemporaryDirectory() as tmp:
        # A compressed CONTCAR, dispatched by basename after the .bz2 suffix is stripped.
        path = Path(tmp) / "CONTCAR.bz2"
        path.write_bytes(bz2.compress(CONTCAR_TEXT.encode("utf-8")))
        data = httk.core.load(str(path))

    print("Format:", data["format"])
    print("Comment:", data["comment"])
    print("Scale:", data["scale"], "Volume:", data["volume"])
    print("Symbols:", data["symbols"], "Counts:", data["counts"])
    print("Cartesian:", data["cartesian"])
    print("Coordinates:", data["coords"])


if __name__ == "__main__":
    main()
