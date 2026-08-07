# Reading VASP POSCAR and CONTCAR files

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

```{literalinclude} ../../examples/parse_poscar.py
:language: python
:lines: 65-
```
