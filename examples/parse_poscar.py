"""Reading VASP POSCAR and CONTCAR files

A POSCAR (or its relaxed sibling CONTCAR) is VASP's structure file: a comment
line, a scaling factor, three lattice vectors, the species and how many atoms of
each, a coordinate-mode line, and then the coordinates. It is a small format
with a surprising number of variants, and `httk.io.read_poscar` handles them
while making one strong promise: **it never converts a number**. Every numeric
field comes back as the verbatim string found in the file, so no rounding
happens at the I/O layer. Turning that mapping into an exact `UnitcellStructure`
is the job of `httk.core.load`.

`read_poscar(source)` accepts a filename (`str` or `os.PathLike`), an open text
stream, or a plain iterable of lines, and returns a neutral, JSON-able mapping:

`format`
    Always `"vasp-poscar"` — a tag identifying which reader produced the mapping.
`comment`
    Line 1, stripped. VASP itself attaches no meaning to it.
`scale` / `volume`
    Line 2, and **exactly one of the two is a string while the other is `None`.**
    A positive universal scaling factor is a multiplier on the lattice vectors
    and lands in `scale`. A *negative* one is VASP's shorthand for "scale the
    cell so its volume equals |this|", which is a completely different physical
    instruction — so it lands in `volume` instead, sign stripped. Collapsing the
    two into one signed number would force every consumer to re-derive the
    distinction, so the reader keeps them apart.
`cell`
    Three rows of three coordinate strings.
`symbols`
    The species line, as `list[str]` — or `None` for a VASP-4 file, where line 6
    holds the atom counts and no species names are recorded anywhere. `None`
    means "this file does not say", not "no species".
`counts`
    Atoms per species, as `list[int]`. These *are* parsed as integers: they are
    counts, not measurements.
`cartesian`
    `True` when the coordinate-mode line begins with C/c/K/k, `False` for
    Direct/fractional. VASP looks only at the first character, and so does this.
`coords`
    One row of three coordinate strings per atom, in file order — `sum(counts)`
    rows.
`selective_dynamics`
    The per-atom T/F flags as booleans when the file declares selective
    dynamics, otherwise `None`.

Two kinds of trailing content are read past and dropped: per-line species labels
after the three coordinates on a coordinate row (some tools write them, VASP
ignores them), and the velocity or predictor-corrector blocks that follow the
coordinates in a CONTCAR. Only `sum(counts)` coordinate rows are ever consumed.

Malformed input raises `ValueError` naming the offending line number, so a
truncated file points at where it stopped making sense rather than failing
somewhere far away.

The last section shows the reader reached through `httk.core.load` instead of
directly. Importing `httk.core` discovers the `httk.registry.io.io` registry
package, which registers this loader under the extensions `.poscar` / `.vasp`
*and* under the exact basenames `POSCAR` / `CONTCAR` — because the canonical
VASP files have no extension at all. `load` therefore dispatches a plain
`CONTCAR`, and also a `CONTCAR.bz2`: the compression suffix is stripped before
the basename is matched, and the file is decompressed transparently on read.

Run this file to see each variant parsed and printed.
"""

import bz2
import tempfile
from pathlib import Path
from typing import Any

import httk.core

from httk.io import read_poscar

# VASP-5: a species line, selective dynamics, and a trailing velocity block.
VASP5_SELECTIVE = """SiO2 slab, top layer relaxed
1.0
4.9134000000000000 0.0000000000000000 0.0000000000000000
0.0000000000000000 4.9134000000000000 0.0000000000000000
0.0000000000000000 0.0000000000000000 5.4052000000000000
Si O
1 2
Selective dynamics
Direct
0.0000000000000000 0.0000000000000000 0.0000000000000000 T T F
0.4130000000000000 0.4130000000000000 0.0000000000000000 F F F
0.5870000000000000 0.5870000000000000 0.5000000000000000 T T T

0.0000000000000000 0.0000000000000000 0.0000000000000000
0.0100000000000000 0.0000000000000000 0.0000000000000000
0.0000000000000000 0.0100000000000000 0.0000000000000000
"""

# VASP-4: no species line at all, cartesian coordinates, and trailing per-line
# labels of the kind some conversion tools emit.
VASP4_CARTESIAN = """Ni3Al, written by a VASP-4 era tool
2.0
1.0000000000 0.0000000000 0.0000000000
0.0000000000 1.0000000000 0.0000000000
0.0000000000 0.0000000000 1.0000000000
3 1
Cartesian
0.0000000000 0.8925000000 0.8925000000 Ni
0.8925000000 0.0000000000 0.8925000000 Ni
0.8925000000 0.8925000000 0.0000000000 Ni
0.0000000000 0.0000000000 0.0000000000 Al
"""

# A negative scaling factor: the number is the target cell VOLUME, not a multiplier.
VOLUME_SCALED = """Cubic He, cell scaled to a fixed volume
-16.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
He
1
Direct
0.0 0.0 0.0
"""


def show_mapping(data: dict[str, Any], *, indent: str = "   ") -> None:
    """Print the neutral mapping, one key per line."""
    for key in ("format", "comment", "scale", "volume", "symbols", "counts", "cartesian"):
        print(f"{indent}{key:<20} {data[key]!r}")
    print(f"{indent}{'cell':<20} (3 lattice-vector rows)")
    for row in data["cell"]:
        print(f"{indent}    {row!r}")
    print(f"{indent}{'coords':<20} ({len(data['coords'])} rows)")
    for row in data["coords"]:
        print(f"{indent}    {row!r}")
    print(f"{indent}{'selective_dynamics':<20} {data['selective_dynamics']!r}")


def show_vasp5_selective() -> None:
    """VASP-5 with selective dynamics: flags read, velocity block ignored."""
    print("== VASP-5: species line, selective dynamics, trailing velocities ==")
    data = read_poscar(VASP5_SELECTIVE.splitlines(keepends=True))
    show_mapping(data)

    print()
    print("The file has", len(VASP5_SELECTIVE.splitlines()), "lines: the coordinates are followed by a blank")
    print("line and a three-row velocity block. Only sum(counts) =", sum(data["counts"]), "coordinate rows are")
    print("ever consumed, so everything after them is read past and dropped.")
    # symbols/counts describe *species*; coords and the flags are *per atom*, in
    # file order. Expanding symbols by counts is how the two line up.
    per_atom = [symbol for symbol, count in zip(data["symbols"], data["counts"]) for _ in range(count)]
    flags = data["selective_dynamics"]
    assert flags is not None
    print()
    for index, (symbol, row) in enumerate(zip(per_atom, flags)):
        free = ", ".join(axis for axis, allowed in zip("xyz", row) if allowed) or "nothing (fully fixed)"
        print(f"   atom {index} ({symbol}): free to move along {free}")
    print("   Coordinates are the verbatim strings from the file: no float() anywhere.")
    print()


def show_vasp4_cartesian() -> None:
    """VASP-4 with cartesian coordinates and trailing labels."""
    print("== VASP-4: no species line, cartesian coordinates, trailing labels ==")
    data = read_poscar(VASP4_CARTESIAN.splitlines(keepends=True))
    show_mapping(data)

    print()
    print("symbols is", data["symbols"], "-- a VASP-4 file records no species names.")
    print("counts still says", data["counts"], "atoms, so the composition's *shape* is known.")
    print("cartesian is", data["cartesian"], "-- the mode line reads 'Cartesian'.")
    print("Each coordinate row in the file ends with a species label ('Ni', 'Al');")
    print("VASP ignores those and so does the reader: only the first three tokens are kept.")
    print()


def show_negative_scale() -> None:
    """A negative scaling factor states a target volume, not a multiplier."""
    print("== Negative scale: the number is a target VOLUME ==")
    data = read_poscar(VOLUME_SCALED.splitlines(keepends=True))
    print("   line 2 of the file:", VOLUME_SCALED.splitlines()[1].strip())
    print("   scale  ->", data["scale"], " (no multiplier was given)")
    print("   volume ->", repr(data["volume"]), "(the sign is stripped; the cell is scaled to this volume)")

    positive = read_poscar(VASP5_SELECTIVE.splitlines(keepends=True))
    print("   for comparison, the VASP-5 file above:")
    print("      scale  ->", repr(positive["scale"]), " volume ->", positive["volume"])
    print("   Exactly one of the two is ever a string; the other is always None.")
    print()


def show_error_reporting() -> None:
    """A malformed file names the line where it stopped making sense."""
    print("== Error reporting ==")
    truncated = "a comment\n1.0\n1 0 0\n0 1 0\n"  # cut off before the third lattice row
    try:
        read_poscar(truncated.splitlines(keepends=True))
    except ValueError as exc:
        print("   truncated file ->", f"{type(exc).__name__}: {exc}")

    bad_flag = VASP5_SELECTIVE.replace("T T F", "T T X")
    try:
        read_poscar(bad_flag.splitlines(keepends=True))
    except ValueError as exc:
        print("   bad T/F flag   ->", f"{type(exc).__name__}: {exc}")
    print()


def show_load_dispatch(directory: Path) -> None:
    """httk.core.load dispatches POSCAR/CONTCAR by basename, compressed or not."""
    print("== httk.core.load: dispatch by extension and by basename ==")
    print("   registered extensions:", httk.core.register.known_extensions())
    print("   registered filenames: ", httk.core.register.known_filenames())
    print("   The canonical VASP files have no extension, hence the basename entries.")

    plain = directory / "CONTCAR"
    plain.write_text(VASP5_SELECTIVE, encoding="utf-8")
    from_plain = httk.core.load(str(plain), raw=True)
    print(f"\n   load({plain.name!r}) -> matched by basename, no extension involved")
    print("      comment:", from_plain["comment"])

    compressed = directory / "CONTCAR.bz2"
    compressed.write_bytes(bz2.compress(VASP5_SELECTIVE.encode("utf-8")))
    from_compressed = httk.core.load(str(compressed), raw=True)
    print(f"   load({compressed.name!r}) -> '.bz2' stripped, then matched by basename")
    print("      symbols:", from_compressed["symbols"], " counts:", from_compressed["counts"])
    print("      decompressed transparently; identical mapping:", from_compressed == from_plain)


def main() -> None:
    show_vasp5_selective()
    show_vasp4_cartesian()
    show_negative_scale()
    show_error_reporting()
    with tempfile.TemporaryDirectory() as tmp:
        show_load_dispatch(Path(tmp))


if __name__ == "__main__":
    main()
