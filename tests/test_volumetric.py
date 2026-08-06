"""Tests for VASP/VESTA volumetric output."""

import io
import pathlib

import pytest

from httk.io.vasp import write_vasp_volumetric

numpy = pytest.importorskip("numpy")


POSCAR_PAYLOAD = {
    "format": "vasp-poscar",
    "comment": "grid",
    "scale": "1.0",
    "volume": None,
    "cell": [["1.0", "0.0", "0.0"], ["0.0", "1.0", "0.0"], ["0.0", "0.0", "1.0"]],
    "symbols": ["Si"],
    "counts": [1],
    "cartesian": False,
    "coords": [["0.0", "0.0", "0.0"]],
    "selective_dynamics": None,
}


def test_volumetric_fortran_order_and_columns(tmp_path: pathlib.Path) -> None:
    grid = numpy.array([[[1.0], [3.0]], [[2.0], [4.0]]])
    path = tmp_path / "density.vasp"
    write_vasp_volumetric(path, POSCAR_PAYLOAD, grid, cols=3)
    expected = (
        "grid\n"
        "1.0\n"
        "1.0 0.0 0.0\n"
        "0.0 1.0 0.0\n"
        "0.0 0.0 1.0\n"
        "Si\n"
        "1\n"
        "Direct\n"
        "0.0 0.0 0.0\n"
        "\n"
        "2 2 1\n"
        "  1.00000000E+00   2.00000000E+00   3.00000000E+00\n"
        "  4.00000000E+00\n"
    )
    assert path.read_text(encoding="utf-8") == expected


def test_volumetric_stream_and_complex_refusal() -> None:
    stream = io.StringIO()
    grid = numpy.arange(2, dtype=numpy.float64).reshape(1, 1, 2)
    write_vasp_volumetric(stream, POSCAR_PAYLOAD, grid)
    assert stream.getvalue().endswith("  0.00000000E+00   1.00000000E+00\n")
    with pytest.raises(ValueError, match=r"\.real.*\.imag"):
        write_vasp_volumetric(io.StringIO(), POSCAR_PAYLOAD, grid.astype(numpy.complex128))
