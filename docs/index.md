# *httk-io*

This site documents specifically the *httk-io* module. For the full
documentation of *httk₂* as a whole, see [docs.httk.org](https://docs.httk.org).

*httk-io* is a *httk₂* module providing file input/output under the namespace
`httk.io`, plus the registry package `httk.registry.io.io` that registers its readers
with *httk-core*. It currently ships a CIF/mCIF parser, reader and writer stack
and a string-preserving VASP POSCAR/CONTCAR reader.

```{admonition} Quick links
:class: tip

- **API reference**: {doc}`reference/index`
- **Reading POSCAR/CONTCAR files**: {doc}`poscar`
- **Reading and writing WAVECAR files**: {doc}`wavecar`
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
block = payload["blocks"][0]
print(payload["header"])
print(block["cell_length_a"])
```

```{toctree}
:maxdepth: 2
:caption: Documentation

poscar
wavecar
examples/index
reference/index
```
