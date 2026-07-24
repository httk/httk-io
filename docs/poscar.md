# Reading VASP POSCAR / CONTCAR files

*httk-io* ships a string-preserving reader for VASP POSCAR/CONTCAR files,
`httk.io.read_poscar`. It parses the file into a neutral, JSON-able mapping whose
numeric fields are kept as the **verbatim strings** found in the file — no
floating-point rounding happens at the I/O layer. Turning that mapping into an
exact `Structure` is the separate job of `httk.atomistic.structure_from_poscar`
(see the *httk-atomistic* docs).

## The neutral mapping

`read_poscar` accepts a filename (`str` or `os.PathLike`), an open text stream,
or a plain iterable of lines. Here we parse a small VASP-5 cell given as a list
of lines:

```python
from httk.io import read_poscar

poscar_lines = [
    "Si primitive\n",
    "1.0\n",
    "0.0 2.7 2.7\n",
    "2.7 0.0 2.7\n",
    "2.7 2.7 0.0\n",
    "Si\n",
    "2\n",
    "Direct\n",
    "0.00 0.00 0.00\n",
    "0.25 0.25 0.25\n",
]

data = read_poscar(poscar_lines)
assert data["format"] == "vasp-poscar"
assert data["comment"] == "Si primitive"
assert data["scale"] == "1.0" and data["volume"] is None
assert data["symbols"] == ["Si"]
assert data["counts"] == [2]
assert data["cartesian"] is False
# Coordinates are preserved exactly as written:
assert data["coords"] == [["0.00", "0.00", "0.00"], ["0.25", "0.25", "0.25"]]
assert data["selective_dynamics"] is None
```

The `scale` and `volume` fields are mutually exclusive: a negative universal
scaling factor on line 2 means its absolute value is the target cell **volume**,
so `scale` is `None` and `volume` carries the (sign-stripped) string. A VASP-4
file (atom counts with no species line) sets `symbols` to `None`. When the file
declares *selective dynamics*, the per-atom `T`/`F` flags are collected into
`selective_dynamics` as booleans; otherwise that field is `None`. Trailing
per-line species labels and velocity blocks are ignored.

## Loader registration and `load`

Importing `httk.core` discovers `httk.handlers.io`, which registers the POSCAR
loader under the extensions `.poscar` / `.vasp` and the exact basenames `POSCAR`
/ `CONTCAR`. `httk.core.load` therefore dispatches these files automatically,
including compressed ones such as `CONTCAR.bz2` (the compression suffix is
stripped to recognize the basename, and the file is decompressed transparently
on read):

```python
import bz2
import tempfile
from pathlib import Path

import httk.core

assert "contcar" in httk.core.register.known_filenames()
assert ".poscar" in httk.core.register.known_extensions()

contcar_text = "He\n1.0\n1 0 0\n0 1 0\n0 0 1\nHe\n1\nDirect\n0 0 0\n"
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "CONTCAR.bz2"
    path.write_bytes(bz2.compress(contcar_text.encode("utf-8")))
    data = httk.core.load(str(path))

assert data["format"] == "vasp-poscar"
assert data["symbols"] == ["He"]
```
