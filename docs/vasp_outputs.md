# Reading VASP output files

The VASP readers keep numeric lexemes as strings. The I/O layer does not parse
numbers into floats; counts and iteration indices are the small structural
integer exceptions.

## POSCAR and CONTCAR

`read_poscar` accepts a path, text stream, or line iterable and returns a
neutral `vasp-poscar` mapping. Path and stream reads include `raw`, the original
decompressed text, when it is available; line iterables set `raw` to `None`.
Saving a payload with `raw` writes that text verbatim, so it takes precedence
over edited tokens and preserves CRLF as well as other formatting:

```python
import httk.core
from httk.io.vasp import read_poscar

text = "c\r\n1.0\r\n1 0 0\r\n0 1 0\r\n0 0 1\r\nHe\r\n1\r\nDirect\r\n0 0 0\r\n"
payload = read_poscar(text.splitlines(keepends=True))
assert payload["raw"] is None

payload = read_poscar(text.splitlines(keepends=True))
assert payload["coords"] == [["0", "0", "0"]]
```

For byte-exact save-back, load a filename (including a compressed filename)
with `raw=True`, then save the returned mapping:

```python
from pathlib import Path

source = Path("POSCAR")
source.write_bytes(text.encode())
payload = httk.core.load(str(source), raw=True)
payload["coords"][0] = ["9", "9", "9"]  # raw still wins
httk.core.save(payload, "POSCAR.copy")
assert Path("POSCAR.copy").read_bytes() == source.read_bytes()
```

## OUTCAR

`read_outcar` returns an `OutcarFile`. Construction checks existence only.
It accepts compressed paths and reopens a fresh forward text stream for each
scan; `close()` is therefore a lifecycle symmetry and owns no persistent handle.
Its `path` property returns the source filename string.
The first access to version, parameters, XC, or POTCAR titles performs a
bounded prologue scan. The first access to finals, completion, issues, frames,
stresses, or elastic moduli performs one full pass and caches only bounded
results. The prologue also exposes `ions_per_type` when VASP prints it.

`frames()` streams complete ionic frames. Each `OutcarFrame` retains cell,
positions, forces, stress, energies, and MD temperature lexemes. `nframes` and
`last_frame` are cached by the full pass; `frame(i)` re-streams the file and is
O(file). `stresses()` returns every six-token `in kB` row. Elastic tables are
available as `ElasticModuliBlock` values for the `TOTAL ELASTIC MODULI`,
`SYMMETRIZED ELASTIC MODULI`, and ionic-relaxation headings.

VASP writes stress as `[xx, yy, zz, xy, yz, zx]` in kBar with compression
positive. `stress_gpa_voigt()` multiplies by `0.1`, flips the sign to make
tension positive, and returns `[xx, yy, zz, yz, xz, xy]`:

```python
from httk.io.vasp import OutcarFrame

frame = OutcarFrame(0, None, None, None, ("1", "2", "3", "4", "5", "6"), None, None, None, None)
assert frame.stress_gpa_voigt() == (-0.1, -0.2, -0.30000000000000004, -0.5, -0.6000000000000001, -0.4)
```

## XDATCAR

`read_xdatcar` returns a lazy `XdatcarFile`. Its bounded header exposes the
comment, scale, cell, symbols, and integer counts. `frames()` supports the
fixed-cell `Direct configuration=` layout and repeated-header variable-cell
layout, yielding coordinate lexemes without retaining all frames. Compressed
paths and context-manager use are supported; an incomplete final block is
dropped and reported in `issues`. The header has `cartesian`, and every frame
mapping has its own `cartesian` boolean from the configuration marker. Repeated-
header frames also expose their validated per-frame `scale`; fixed-cell frames
use `None`. `XdatcarFile.path` returns the source filename string.

## OSZICAR and POTCAR summary

`read_oszicar` groups electronic iterations with the following ionic summary.
Electronic values and ionic energies remain strings. A trailing electronic
block becomes an entry whose ionic fields are `None`. MD summaries may omit
`dE`; that missing lexeme is retained as `None` rather than inferred.

`read_potcar_summary` extracts only header metadata (`TITEL`, `ZVAL`, `POMASS`,
`ENMAX`, and `LEXCH`) into one mapping per potential. It supports concatenated
`POTCAR.summary` headers and never retains or exposes the full POTCAR text.

## VASPOutputs

`VASPOutputs(directory)` lazily probes `POSCAR`, `CONTCAR`, `OUTCAR`,
`XDATCAR`, `OSZICAR`, and `POTCAR`/`POTCAR.summary`, including every compression
suffix registered by *httk-core*. Missing files return `None`. Payload-returning
properties hold no handles; the composite owns and closes the `OutcarFile` and
`XdatcarFile` objects it constructs:

```python
from httk.io.vasp import VASPOutputs

with VASPOutputs("calculation") as outputs:
    assert outputs.poscar is not None
    assert outputs.outcar is None or outputs.outcar.version_string.startswith("vasp.")
```
