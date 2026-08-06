"""Tests for the lazy XDATCAR reader."""

import bz2
from pathlib import Path

import pytest

from httk.io.vasp import XdatcarFile, read_xdatcar

HEADER = """Synthetic XDATCAR comment with spaces
  1.000000
  1.000000  0.000000  0.000000
  0.000000  2.000000  0.000000
  0.000000  0.000000  3.000000
  Si  O
  1  2
"""

FIXED = (
    HEADER
    + """Direct configuration=     1
  0.100000  0.200000  0.300000
  0.400000  0.500000  0.600000
  0.700000  0.800000  0.900000
Direct configuration=     2
  0.110000  0.210000  0.310000
  0.410000  0.510000  0.610000
  0.710000  0.810000  0.910000
"""
)


def variable_header(cell: str, scale: str = "1.000000") -> str:
    return HEADER.splitlines(keepends=True)[0] + f"{scale}\n" + cell + "\n  Si  O\n  1  2\n"


CELL_1 = "  4.000000  0.000000  0.000000\n  0.000000  4.000000  0.000000\n  0.000000  0.000000  4.000000"
CELL_2 = "  5.000000  0.000000  0.000000\n  0.000000  5.000000  0.000000\n  0.000000  0.000000  5.000000"
NPT = (
    HEADER
    + variable_header(CELL_1)
    + "Direct configuration= 1\n"
    + """0 0 0
0 0 0
0 0 0
"""
    + variable_header(CELL_2, "2.000000")
    + "Direct configuration= 2\n"
    + """1 1 1
1 1 1
1 1 1
"""
)


def test_fixed_cell_header_and_frames(tmp_path: Path) -> None:
    path = tmp_path / "XDATCAR"
    path.write_text(FIXED, encoding="utf-8")
    xdatcar = XdatcarFile(path)
    assert xdatcar.path == str(path)
    assert xdatcar.comment == "Synthetic XDATCAR comment with spaces"
    assert xdatcar.scale == "1.000000"
    assert xdatcar.cell == (
        ("1.000000", "0.000000", "0.000000"),
        ("0.000000", "2.000000", "0.000000"),
        ("0.000000", "0.000000", "3.000000"),
    )
    assert xdatcar.symbols == ("Si", "O")
    assert xdatcar.counts == (1, 2)
    frames = tuple(xdatcar.frames())
    assert frames[0] == {
        "index": 0,
        "cell": None,
        "cartesian": False,
        "scale": None,
        "coords": (
            ("0.100000", "0.200000", "0.300000"),
            ("0.400000", "0.500000", "0.600000"),
            ("0.700000", "0.800000", "0.900000"),
        ),
    }
    assert frames[1]["index"] == 1
    assert xdatcar.nframes == 2
    assert xdatcar.issues == ()


def test_cartesian_mode_is_exposed(tmp_path: Path) -> None:
    path = tmp_path / "XDATCAR"
    path.write_text(HEADER + "Cartesian configuration= 1\n0 0 0\n0 0 0\n0 0 0\n", encoding="utf-8")
    xdatcar = XdatcarFile(path)
    assert xdatcar.cartesian is True
    assert next(xdatcar.frames())["cartesian"] is True


def test_npt_repeated_headers_and_compression(tmp_path: Path) -> None:
    path = tmp_path / "XDATCAR.bz2"
    path.write_bytes(bz2.compress(NPT.encode("utf-8")))
    xdatcar = XdatcarFile(path)
    frames = tuple(xdatcar.frames())
    assert frames[0]["cell"] == (
        ("4.000000", "0.000000", "0.000000"),
        ("0.000000", "4.000000", "0.000000"),
        ("0.000000", "0.000000", "4.000000"),
    )
    assert frames[1]["cell"] == (
        ("5.000000", "0.000000", "0.000000"),
        ("0.000000", "5.000000", "0.000000"),
        ("0.000000", "0.000000", "5.000000"),
    )
    assert frames[0]["scale"] == "1.000000"
    assert frames[1]["scale"] == "2.000000"
    assert xdatcar.nframes == 2


@pytest.mark.parametrize(
    ("name", "text", "message"),
    (
        (
            "bad_scale",
            NPT.replace("2.000000\n  5.000000", "bad-scale\n  5.000000", 1),
            "malformed XDATCAR",
        ),
        (
            "bad_symbols",
            NPT.replace("  Si  O\n  1  2\nDirect configuration= 2", "  Si  C\n  1  2\nDirect configuration= 2"),
            "symbols/counts differ",
        ),
        (
            "missing_marker",
            NPT.replace("Direct configuration= 2", "not a configuration marker"),
            "not followed by configuration marker",
        ),
        (
            "discontinuous",
            FIXED.replace("configuration=     2", "configuration=     4"),
            "discontinuous configuration index",
        ),
    ),
)
def test_xdatcar_rejects_bad_repeated_headers_or_indices(tmp_path: Path, name: str, text: str, message: str) -> None:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    xdatcar = XdatcarFile(path)
    assert xdatcar.path == str(path)
    tuple(xdatcar.frames())
    assert any(message in issue for issue in xdatcar.issues)


def test_truncated_final_block_is_dropped_and_reported(tmp_path: Path) -> None:
    path = tmp_path / "XDATCAR"
    path.write_text(FIXED + "Direct configuration= 3\n0 0 0\n", encoding="utf-8")
    xdatcar = XdatcarFile(path)
    assert len(tuple(xdatcar.frames())) == 2
    assert len(xdatcar.issues) == 1
    assert "truncated final coordinate block" in xdatcar.issues[0]
    assert xdatcar.nframes == 2


def test_reader_lifecycle_and_registry(tmp_path: Path) -> None:
    path = tmp_path / "XDATCAR"
    path.write_text(FIXED, encoding="utf-8")
    payload = read_xdatcar(path)
    assert payload["format"] == "vasp-xdatcar"
    with XdatcarFile(path) as xdatcar:
        assert not xdatcar.closed
    assert xdatcar.closed
    with pytest.raises(ValueError, match="closed XDATCAR"):
        _ = xdatcar.nframes
