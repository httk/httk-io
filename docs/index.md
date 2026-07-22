# *httk-io*

*httk-io* is a *httk v2* module providing file input/output under the namespace
`httk.io`, plus the handler package `httk.handlers.io` that registers its loaders
with *httk-core*. It currently ships a CIF/mCIF parser, reader and writer stack.

```{admonition} Quick links
:class: tip

- **API reference**: {doc}`reference/index`
````

## Install

Preferably work in a Python virtual environment, then do:
```bash
git clone https://github.com/httk/httk-io
cd httk-io
python -m pip install -e .
```

## Usage example

Importing `httk.core` discovers `httk.handlers.io` and registers the `.cif`
loader, so `httk.core.load` can dispatch a CIF file to *httk-io*. The loader
returns the raw parsed CIF data as a `(data_blocks, header)` tuple:

```python
import httk.core  # discovery registers the ".cif" loader

data_blocks, header = httk.core.load("structure.cif")
name, block = data_blocks[0]
print(block["cell_length_a"])
```

```{toctree}
:maxdepth: 2
:caption: Documentation

reference/index
```
