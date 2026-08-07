# Reading, inspecting and writing CIF files

The Crystallographic Information File (CIF) format is the standard exchange
format for crystal structures. *httk-io* provides a low-level CIF stack that
treats a CIF as what it literally is — a sequence of named data blocks, each a
set of tag/value pairs plus tabular `loop_` sections — and refuses to guess what
any of it *means*. Interpreting the tags (turning a cell and an asymmetric unit
into a structure) is a separate, higher-level job.

Three pieces make up that stack, and this example walks through all of them.

## `read_cif(source)`

Parses a CIF and returns a `(data_blocks, header)` tuple. `data_blocks` is a
list of `(name, block)` pairs — CIF permits several `data_` blocks in one file,
so this is a list, not a dict. `header` is the run of comment lines at the very
top of the file, kept verbatim so that it survives a round trip.

A `block` is a flat mapping from *lowercased tag name without its leading
underscore* to value. Two rules make that mapping worth knowing:

- **Values are strings.** `_cell_length_a 5.6402(3)` becomes the string
  `"5.6402(3)"`, not a float. Nothing is rounded, nothing is discarded, and the
  parenthesised standard uncertainty is still there for whoever wants it.
- **Loops become columns.** Each `loop_` contributes one list-valued entry per
  column, all of equal length, plus a bookkeeping key `loop_0`, `loop_1`, …
  listing that loop's column names in file order. Those bookkeeping keys are
  what let `write_cif` reconstruct the loops, and what let a reader tell a
  one-row loop from a plain tag.

`source` may be a filename (`str` or `os.PathLike`), an already-open text
stream, or any iterable of lines. Filenames are opened through
`httk.core.TextstreamFileView`, so `structure.cif.bz2` and `structure.cif.gz`
are decompressed transparently.

## `parse_cif_float(token)` and `parse_cif_int(token)`

These convert one of the preserved strings into a number, honouring the numeric
conventions CIF actually uses:

- `"5.6402(3)"` — a value with a standard uncertainty (esd) on its last digits.
  The plain call returns just the central value `5.6402`; passing `meta=True`
  returns `(value, {"esd": ..., "precision": ...})`, where `esd` is the
  uncertainty in the value's own units (`3/10000` here) and `precision` is the
  precision implied by how the number was *written* (`0.0001` for four
  decimals). The two are independent: a value with no esd still has a stated
  precision.
- `"?"` — "unknown", a legal CIF value. It parses to `None` rather than raising,
  so a missing measurement stays representable. Its sibling `"."` ("not
  applicable") is an error instead: there is no number to return.
- `"1/3"` — a fraction, common in symmetry operations. Parsed exactly and
  reported with precision `None`, i.e. a stated value rather than a rounded decimal.

`parse_cif_int` follows the same conventions and returns the integral central
value: `"295(2)"` is `295`, `"3E2"` is `300`. It is strict by default, so a
genuinely fractional token such as `"1.5"` raises `ValueError` instead of
silently rounding.

## `write_cif(target, data_blocks, header)`

Writes the parsed structure back out: the header comments first, then each block
as `data_<name>` followed by its tags and its reconstructed loops. Values the
writer cannot recognise as bare CIF numbers are quoted, which is why an
esd-bearing string such as `5.6402(3)` comes back out as `'5.6402(3)'` —
re-reading strips the quotes, so the value is unchanged. Demonstrating that is
the point of the round-trip section below: read → write → read yields an
identical mapping.

The final section shows the same file arriving through `httk.core.load`, which
dispatches `.cif` (and `.cif.bz2`) to *httk-io* because importing `httk.core`
discovers the `httk.registry.io.io` registry package. By default, `load(path)` returns
the native `ASUStructure` when *httk-atomistic* is installed. This example uses
`load(path, raw=True)`, which returns the neutral payload mapping instead. Blocks
that are not structures are reported rather than raising, since a CIF may hold
anything.

Run this file to see every step printed.

```{literalinclude} ../../examples/parse_cif.py
:language: python
:lines: 79-
```
