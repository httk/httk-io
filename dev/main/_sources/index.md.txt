# *httk-io*

This site documents specifically the *httk-io* module. For the full
documentation of *httk₂* as a whole, see [docs.httk.org](https://docs.httk.org).

*httk-io* is a *httk₂* module providing file input/output under the namespace
`httk.io`, plus the registry package `httk.registry.io.io` that registers its readers
with *httk-core*. It currently ships a CIF/mCIF parser, reader and writer stack;
string-preserving VASP POSCAR/CONTCAR readers plus the OUTCAR/XDATCAR/OSZICAR/POTCAR
output readers; a numpy-backed WAVECAR reader/writer and volumetric-grid writer; and
the OPTIMADE trajectory JSON Lines holding format.

```{admonition} Quick links
:class: tip

- **API reference**: {doc}`reference/index`
- **Reading POSCAR/CONTCAR files**: {doc}`poscar`
- **Reading VASP output files**: {doc}`vasp_outputs`
- **Reading and writing WAVECAR files**: {doc}`wavecar`
- **Trajectory JSON Lines**: {doc}`trajectory_jsonl`
- **Runnable examples**: {doc}`examples/index`
````

## Install

Preferably work in a Python virtual environment, then do:
```bash
git clone https://github.com/httk/httk-io
cd httk-io
python -m pip install -e .
```

## Usage example

Importing `httk.core` discovers `httk.registry.io.io` and registers the `.cif`
loader, so `httk.core.load` can dispatch a CIF file to *httk-io*. The loader
returns a neutral parsed CIF payload when called with `raw=True`:

```python
import httk.core  # discovery registers the ".cif" loader

payload = httk.core.load("structure.cif", raw=True)
block = payload["blocks"][0]           # one neutral asymmetric-unit mapping per structural block
print(payload["header"])               # the file's leading comment lines
print(block["symbols"])                # e.g. ["Na", "Cl"]
print(block["cell_parameters_exact"])  # ("a", "b", "c", "alpha", "beta", "gamma") as verbatim tokens
```

## CIF loading

The neutral CIF payload is a mapping with `format` `"cif"` (`"mcif"` for magnetic
CIFs), a `blocks` list holding one asymmetric-unit mapping per structural block,
`unparsed` reasons for blocks that have atom sites but cannot be interpreted, and
the verbatim `header`. Numeric values are kept as strings; where the file carries
`_httk_*_exact` companion tags those exact tokens are preferred, so no precision is
lost at the I/O layer. `.cif`, `.mcif` and their compressed forms (`.cif.gz`,
`.cif.bz2`) all dispatch here.

Two conveniences smooth over real-world files:

- **Inferred element symbols.** `_atom_site_type_symbol` is optional in the CIF
  core dictionary. When it is absent, each site's element is inferred from the
  leading element run of its `_atom_site_label` (`"MgM1"` → `Mg`), with a
  `RuntimeWarning`. A label whose prefix names no element is not guessed at: its
  block cannot be interpreted, so `load` omits it from `blocks` and records the
  reason in `unparsed` (the underlying parser raises a `ValueError` that the
  loader catches per block).
- **Autocorrect.** Passing `autocorrect=True` to `load` (or to `read_cif` /
  `read_cif_asus`) drops a malformed *auxiliary* loop — one whose column counts do
  not line up and whose tags are not a protected structural family — warning about
  each repair instead of refusing the file, and stamps `autocorrect=True` on the
  payload. Without it, such a loop is a hard `ValueError`.

The runnable {doc}`examples/parse_cif` walks the lower-level `read_cif`,
`parse_cif_float` / `parse_cif_int` and `write_cif` API in full.

```{toctree}
:maxdepth: 2
:caption: Documentation

poscar
vasp_outputs
wavecar
trajectory_jsonl
examples/index
reference/index
```
