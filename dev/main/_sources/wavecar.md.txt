# Reading and writing VASP WAVECAR files

*httk-io* provides the numpy-backed binary WAVECAR layer through
`httk.io.vasp`. It reads the small metadata headers eagerly and reads one
coefficient vector at a time, which keeps the large coefficient data on disk
until requested.

## Installation and dispatch

WAVECAR support requires the optional numpy extra:

```bash
python -m pip install -e '.[numpy]'
```

Importing `httk.core` discovers the WAVECAR reader and writer. Files named
`WAVECAR` and files with the `.wavecar` extension are registered. The core
loader returns a neutral payload:

```python
import httk.core

payload = httk.core.load("WAVECAR", raw=True)
assert payload["format"] == "vasp-wavecar"
source = payload["wavecar"]
```

The `wavecar` value is a `WavecarSource` contract, not an atomistic-domain
object. It exposes `nspins`, `nkpts`, `nbands`, `encut`, `cell`, `kpoints`,
`eigenvalues`, `occupations`, `nplanewaves`, `double_precision`,
and `coefficients(spin, kpt, band)`. `record_length` is specific to an open
`WavecarFile` and is not part of the common source contract. The atomistic
layer can provide the same contract for in-memory data without depending on
the concrete `WavecarFile` class.

All indices are zero-based:

```python
coefficients = source.coefficients(spin=0, kpt=0, band=0)
```

`WavecarFile` is also a context manager. Its metadata arrays have shapes
`(nkpts, 3)`, `(nspins, nkpts, nbands)`, and `(nkpts,)` for k-points,
eigenvalues/occupations, and plane-wave counts respectively.

## Compression and writing

WAVECAR is a random-access binary format. Compressed paths are deliberately
refused: streaming decompression cannot support seeking to arbitrary
spin/k-point/band coefficient records. Decompress a file on disk before
reading it. Writing likewise requires a binary filesystem path rather than a
compressed or text stream.

```python
from httk.io.vasp import read_wavecar, write_wavecar

payload = read_wavecar("WAVECAR")
write_wavecar("WAVECAR.copy", payload)
```

`write_wavecar` accepts the neutral `vasp-wavecar` payload and any source that
implements the `WavecarSource` contract. It preserves single- versus
double-precision coefficient storage according to `double_precision`.

## VASP/VESTA volumetric output

`write_vasp_volumetric` writes POSCAR content followed by a three-dimensional,
real-valued grid in Fortran order. It is the low-level writer used when a
wavefunction or other scalar field needs to be opened by VESTA; pass
`grid.real` or `grid.imag` explicitly for a complex array.

```python
from httk.io.vasp import write_vasp_volumetric

write_vasp_volumetric("wave_r.vasp", poscar_payload, wave.real)
write_vasp_volumetric("wave_i.vasp", poscar_payload, wave.imag)
```
