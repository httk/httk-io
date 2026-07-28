"""Reading, inspecting and writing CIF files

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
discovers the `httk.handlers.io` handler package. Note that `load` returns the
*interpreted* payload — a `format` tag plus one asymmetric unit per structural
data block — rather than the token tree shown above, so that a CIF can be handed
straight to `httk.atomistic.load_structure`. Blocks that are not structures are
reported rather than raising, since a CIF may hold anything.

Run this file to see every step printed.
"""

import bz2
import tempfile
from pathlib import Path
from typing import Any

import httk.core

from httk.io.cif.cif_parser import parse_cif_float, parse_cif_int
from httk.io.cif.cif_reader import read_cif
from httk.io.cif.cif_writer import write_cif

# A small, deliberately messy CIF: two header comment lines, mixed-case tags,
# quoted strings, values carrying standard uncertainties, an unknown value ('?'),
# and two loops of different widths.
NACL_CIF = """#  Sodium chloride, rock-salt structure
#  Toy data for the httk-io example -- not from a real refinement.
data_nacl
_chemical_formula_sum             'Na Cl'
_cell_length_a                    5.6402(3)
_cell_length_b                    5.6402(3)
_cell_length_c                    5.6402(3)
_cell_angle_alpha                 90.0
_cell_volume                      179.43(2)
_cell_measurement_temperature     295(2)
_symmetry_space_group_name_H-M    'F m -3 m'
_symmetry_Int_Tables_number       225
_diffrn_ambient_pressure          ?
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Na1 Na 0.0 0.0 0.0 1.0
Cl1 Cl 0.5 0.5 0.5 1.0
loop_
_symmetry_equiv_pos_as_xyz
'x, y, z'
'-x, -y, z'
'x+1/2, y+1/2, z'
"""


def show_structure(data_blocks: list[Any], header: str) -> dict[str, Any]:
    """What `read_cif` hands back: the header, the blocks, the tags and the loops."""
    print("== read_cif: the (data_blocks, header) tuple ==")
    print("header lines:         ", header.splitlines())
    print("number of data blocks:", len(data_blocks))

    name, block = data_blocks[0]
    print("first block name:     ", repr(name), "(the file's 'data_nacl', lowercased)")

    # A block is flat: loop bookkeeping keys name the loop columns, and every
    # remaining key that is not itself a loop column is a plain tag.
    loop_keys = [key for key in block if key.startswith("loop_")]
    loop_columns = {column for key in loop_keys for column in block[key]}
    tags = [key for key in block if key not in loop_keys and key not in loop_columns]

    print("\n-- plain tags (every value is a string) --")
    for tag in tags:
        print(f"   {tag:<32} {block[tag]!r}")

    print("\n-- loops (one list-valued key per column, plus the loop_N key) --")
    for key in loop_keys:
        columns = block[key]
        print(f"   {key} -> {columns}")
        for row in zip(*(block[column] for column in columns)):
            print("        ", "  ".join(repr(value) for value in row))

    print("\nTag names are normalised: the file writes _symmetry_space_group_name_H-M,")
    print("the block stores", repr("symmetry_space_group_name_h-m"), "-- lowercased, underscore stripped.")
    print()
    return block


def show_numeric_fields(block: dict[str, Any]) -> None:
    """CIF numbers carry uncertainties and special values; the parsers know it."""
    print("== parse_cif_float / parse_cif_int: CIF numeric conventions ==")

    a_token = block["cell_length_a"]
    a_value, a_meta = parse_cif_float(a_token, meta=True)
    print(f"cell_length_a    {a_token!r}")
    print(f"    value      = {a_value}")
    print(f"    esd        = {a_meta['esd']}     (the '(3)' applies to the last digit: 0.0003 Angstrom)")
    print(f"    precision  = {a_meta['precision']}     (implied by writing four decimals)")
    print(f"    plain call, no meta: parse_cif_float({a_token!r}) -> {parse_cif_float(a_token)}")

    t_token = block["cell_measurement_temperature"]
    t_value, t_meta = parse_cif_float(t_token, meta=True)
    print(f"temperature      {t_token!r} -> {t_value} +/- {t_meta['esd']} K")

    angle_token = block["cell_angle_alpha"]
    angle_value, angle_meta = parse_cif_float(angle_token, meta=True)
    print(
        f"angle_alpha      {angle_token!r} -> {angle_value}, esd {angle_meta['esd']}, "
        f"precision {angle_meta['precision']} (no esd, but a stated precision)"
    )

    unknown_token = block["diffrn_ambient_pressure"]
    print(f"unknown value    {unknown_token!r} -> {parse_cif_float(unknown_token)}   (CIF '?' means 'unknown')")
    print("'.' (not applicable) has no numeric value at all: ", end="")
    try:
        parse_cif_float(".")
    except Exception as exc:  # the parser raises a bare Exception here; showing it is the point
        print(f"{type(exc).__name__}: {exc}")

    third = parse_cif_float("1/3")
    _, third_meta = parse_cif_float("1/3", meta=True)
    print(f"fraction         '1/3' -> {third} (precision {third_meta['precision']}: exact, not a rounded decimal)")

    print("\n-- parse_cif_int: the integral central value --")
    tables_number = block["symmetry_int_tables_number"]
    print(f"    {tables_number!r:<8} -> {parse_cif_int(tables_number)} (the space-group number, as an int)")
    print(f"    {'295(2)'!r:<8} -> {parse_cif_int('295(2)')} (the esd is dropped)")
    print(f"    {'3E2'!r:<8} -> {parse_cif_int('3E2')} (exponents are honoured)")
    print(f"    {'1.5'!r:<8} -> ", end="")
    try:
        parse_cif_int("1.5")
    except ValueError as exc:
        print(f"ValueError: {exc}")
    print()


def show_roundtrip(directory: Path, data_blocks: list[Any], header: str) -> None:
    """read -> write -> read gives back an identical mapping."""
    print("== write_cif: writing the parsed data back out ==")
    out_path = directory / "nacl-out.cif"
    write_cif(str(out_path), data_blocks, header)

    print(f"-- {out_path.name} --")
    for line in out_path.read_text(encoding="utf-8").splitlines():
        print("   |", line)

    reread_blocks, reread_header = read_cif(str(out_path))
    print()
    print("re-read: identical header:      ", reread_header == header)
    print("re-read: identical block name:  ", reread_blocks[0][0] == data_blocks[0][0])
    print("re-read: identical block mapping:", dict(reread_blocks[0][1]) == dict(data_blocks[0][1]))
    print()
    print("The written file quotes the values the writer cannot see as bare numbers")
    print("('5.6402(3)', '?'); re-reading strips the quotes, so nothing has changed.")
    print()


def show_load_dispatch(directory: Path) -> None:
    """The interpreting reader, reached through httk.core.load — including compressed."""
    print("== httk.core.load: the registered '.cif' loader ==")
    print("registered extensions:", httk.core.register.known_extensions())
    print("Importing httk.core discovers httk.handlers.io, which registers these.")

    compressed = directory / "nacl.cif.bz2"
    compressed.write_bytes(bz2.compress(NACL_CIF.encode("utf-8")))
    payload = httk.core.load(str(compressed), raw=True)

    print(f"load({compressed.name!r}) -> format={payload['format']!r}; .bz2 decompressed transparently")
    print("   header first line:  ", payload["header"].splitlines()[0])
    print("   structural blocks:  ", len(payload["blocks"]))
    for item in payload["unparsed"]:
        print(f"   not a structure:     block {item['block']!r} ({item['reason']})")
    print("This block does not carry everything a structure needs, so no asymmetric unit")
    print("could be built from it. Loading still succeeds and says why, because a CIF is")
    print("not obliged to describe a crystal at all.")

    # The tags themselves stay reachable through the low-level tokenizer.
    data_blocks, _header = read_cif(str(compressed))
    print("   read_cif tags:       cell_length_a =", data_blocks[0][1]["cell_length_a"])


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        cif_path = directory / "nacl.cif"
        cif_path.write_text(NACL_CIF, encoding="utf-8")

        data_blocks, header = read_cif(str(cif_path))

        block = show_structure(data_blocks, header)
        show_numeric_fields(block)
        show_roundtrip(directory, data_blocks, header)
        show_load_dispatch(directory)


if __name__ == "__main__":
    main()
