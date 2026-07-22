# httk-io

*httk-io* is the file input/output module for *httk₂*. It provides:

- the `httk.io` package within the PEP 420 native `httk` namespace, currently a
  CIF/mCIF parser, reader and writer stack; and
- the `httk.handlers.io` package, which registers the `.cif` loader with
  `httk-core`'s plugin registry so files can be read via `httk.core.load`.

Most users should install the [`httk2`](https://github.com/httk/httk2)
metapackage, which selects a useful set of httk modules:

```console
pip install httk2
```

Install only this module (it depends on `httk-core` and `numpy`) with:

```console
pip install httk-io
```

## Usage

Importing `httk.core` discovers `httk.handlers.io` and registers the `.cif`
loader. `httk.core.load` then dispatches on the file extension and returns the
raw parsed CIF data as a `(data_blocks, header)` tuple:

```python
import httk.core  # discovery registers the ".cif" loader

data_blocks, header = httk.core.load("structure.cif")
name, block = data_blocks[0]
print(block["cell_length_a"])
```

Development and release instructions are in the
[`RELEASING.md`](https://github.com/httk/httk-io/blob/main/RELEASING.md) guide.
