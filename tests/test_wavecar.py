"""Tests for the neutral VASP WAVECAR binary layer."""

import importlib
import io
import os
import pathlib
import struct
import subprocess
import sys

import httk.core
import pytest

try:
    import numpy
except ImportError:  # pragma: no cover - exercised by the subprocess test below
    numpy = None

requires_numpy = pytest.mark.skipif(numpy is None, reason="WAVECAR tests require numpy")
if numpy is not None:
    from httk.io.vasp.wavecar import WavecarFile, _write_wavecar_payload, read_wavecar


def _record(values: list[float], record_length: int) -> bytes:
    body = struct.pack("=" + "d" * len(values), *values)
    assert len(body) <= record_length
    return body + b"\0" * (record_length - len(body))


def _golden_data() -> dict:
    return {
        "cell": numpy.array([[3.0, 0.1, 0.2], [0.3, 4.0, 0.4], [0.5, 0.6, 5.0]], dtype=numpy.float64),
        "kpoints": numpy.array([[0.0, 0.0, 0.0], [0.25, 0.5, 0.75]], dtype=numpy.float64),
        "eigenvalues": numpy.array([[[1.1, 2.2], [3.3, 4.4]]], dtype=numpy.float64),
        "occupations": numpy.array([[[1.0, 0.5], [0.25, 0.0]]], dtype=numpy.float64),
        "nplanewaves": numpy.array([3, 4], dtype=numpy.int64),
        "coefficients": {
            (0, 0, 0): numpy.array([1 + 2j, 3 + 4j, 5 + 6j], dtype=numpy.complex64),
            (0, 0, 1): numpy.array([7 + 8j, 9 + 10j, 11 + 12j], dtype=numpy.complex64),
            (0, 1, 0): numpy.array([13 + 14j, 15 + 16j, 17 + 18j, 19 + 20j], dtype=numpy.complex64),
            (0, 1, 1): numpy.array([21 + 22j, 23 + 24j, 25 + 26j, 27 + 28j], dtype=numpy.complex64),
        },
    }


def _golden_bytes(data: dict, record_length: int, dtype) -> bytes:
    records = [_record([record_length, 1, 45200], record_length)]
    records.append(_record([2, 2, 350.0, *data["cell"].reshape(-1)], record_length))
    for kpt in range(2):
        header = [3 if kpt == 0 else 4, *data["kpoints"][kpt]]
        for band in range(2):
            header.extend([data["eigenvalues"][0, kpt, band], 0.0, data["occupations"][0, kpt, band]])
        records.append(_record(header, record_length))
        for band in range(2):
            body = numpy.asarray(data["coefficients"][(0, kpt, band)], dtype=dtype).tobytes()
            records.append(body + b"\0" * (record_length - len(body)))
    return b"".join(records)


class Source:
    def __init__(self, data: dict, *, double_precision: bool = False, nbands: int = 2) -> None:
        self.nspins = 1
        self.nkpts = 2 if nbands == 2 else 1
        self.nbands = nbands
        self.encut = 350.0
        self.cell = data["cell"]
        self.kpoints = data["kpoints"] if nbands == 2 else numpy.zeros((1, 3))
        self.eigenvalues = (
            data["eigenvalues"] if nbands == 2 else numpy.arange(nbands, dtype=numpy.float64).reshape(1, 1, nbands)
        )
        self.occupations = data["occupations"] if nbands == 2 else numpy.ones((1, 1, nbands), dtype=numpy.float64)
        self.nplanewaves = data["nplanewaves"] if nbands == 2 else numpy.array([1], dtype=numpy.int64)
        self.double_precision = double_precision
        self.record_length = 0
        self._coefficients = (
            data["coefficients"]
            if nbands == 2
            else {(0, 0, band): numpy.array([band + 1j * band], dtype=numpy.complex128) for band in range(nbands)}
        )

    def coefficients(self, spin: int, kpt: int, band: int):
        return self._coefficients[(spin, kpt, band)]


@requires_numpy
def test_golden_read(tmp_path: pathlib.Path) -> None:
    data = _golden_data()
    path = tmp_path / "WAVECAR"
    path.write_bytes(_golden_bytes(data, 128, numpy.complex64))

    with WavecarFile(path) as wavecar:
        assert wavecar.nspins == 1
        assert wavecar.nkpts == 2
        assert wavecar.nbands == 2
        assert wavecar.encut == 350.0
        numpy.testing.assert_array_equal(wavecar.cell, data["cell"])
        numpy.testing.assert_array_equal(wavecar.kpoints, data["kpoints"])
        numpy.testing.assert_array_equal(wavecar.eigenvalues, data["eigenvalues"])
        numpy.testing.assert_array_equal(wavecar.occupations, data["occupations"])
        numpy.testing.assert_array_equal(wavecar.nplanewaves, data["nplanewaves"])
        assert wavecar.double_precision is False
        assert wavecar.record_length == 128
        for key, expected in data["coefficients"].items():
            numpy.testing.assert_array_equal(wavecar.coefficients(*key), expected)


@requires_numpy
def test_gamma_half_hint_round_trips_without_file_semantics(tmp_path: pathlib.Path) -> None:
    data = _golden_data()
    path = tmp_path / "WAVECAR"
    path.write_bytes(_golden_bytes(data, 128, numpy.complex64))
    for gamma_half in (None, "x", "z"):
        payload = read_wavecar(path, gamma_half=gamma_half)
        try:
            assert payload["gamma_half"] == gamma_half
        finally:
            payload["wavecar"].close()
    with pytest.raises(ValueError, match="gamma_half"):
        read_wavecar(path, gamma_half="y")


@requires_numpy
def test_golden_write_is_byte_exact(tmp_path: pathlib.Path) -> None:
    data = _golden_data()
    source = Source(data)
    path = tmp_path / "WAVECAR"
    _write_wavecar_payload(path, {"format": "vasp-wavecar", "wavecar": source, "gamma_half": "z"})

    expected = _golden_bytes(data, 96, numpy.complex64)
    assert path.read_bytes() == expected


@requires_numpy
@pytest.mark.parametrize("double_precision", (False, True))
def test_round_trip_both_precisions(tmp_path: pathlib.Path, double_precision: bool) -> None:
    data = _golden_data()
    if double_precision:
        data["coefficients"] = {key: value.astype(numpy.complex128) for key, value in data["coefficients"].items()}
    source = Source(data, double_precision=double_precision)
    path = tmp_path / "WAVECAR"
    _write_wavecar_payload(path, {"format": "vasp-wavecar", "wavecar": source})
    with WavecarFile(path) as result:
        assert result.double_precision is double_precision
        numpy.testing.assert_array_equal(result.cell, source.cell)
        numpy.testing.assert_array_equal(result.kpoints, source.kpoints)
        numpy.testing.assert_array_equal(result.eigenvalues, source.eigenvalues)
        numpy.testing.assert_array_equal(result.occupations, source.occupations)
        numpy.testing.assert_array_equal(result.nplanewaves, source.nplanewaves)
        for key in source._coefficients:
            numpy.testing.assert_array_equal(result.coefficients(*key), source.coefficients(*key))
    assert result.closed


@requires_numpy
def test_header_alignment_regression(tmp_path: pathlib.Path) -> None:
    data = _golden_data()
    source = Source(data, double_precision=True, nbands=3)
    path = tmp_path / "WAVECAR"
    _write_wavecar_payload(path, {"format": "vasp-wavecar", "wavecar": source})
    with WavecarFile(path) as result:
        assert result.record_length == 112  # 104-byte header rounded up to complex128 alignment.
        for band in range(3):
            numpy.testing.assert_array_equal(result.coefficients(0, 0, band), source.coefficients(0, 0, band))


@requires_numpy
def test_wavecar_validations(tmp_path: pathlib.Path) -> None:
    data = _golden_data()
    raw = bytearray(_golden_bytes(data, 128, numpy.complex64))

    bad_rtag = bytearray(raw)
    struct.pack_into("=d", bad_rtag, 16, 999.0)
    bad_rtag_path = tmp_path / "bad-rtag"
    bad_rtag_path.write_bytes(bad_rtag)
    with pytest.raises(ValueError, match="RTAG.*999"):
        WavecarFile(bad_rtag_path)

    non_integral = bytearray(raw)
    struct.pack_into("=d", non_integral, 128, 1.5)
    non_integral_path = tmp_path / "non-integral"
    non_integral_path.write_bytes(non_integral)
    with pytest.raises(ValueError, match="k-points.*positive integer"):
        WavecarFile(non_integral_path)

    truncated_path = tmp_path / "truncated"
    truncated_path.write_bytes(raw[:-1])
    with pytest.raises(ValueError, match="truncated"):
        WavecarFile(truncated_path)

    path = tmp_path / "WAVECAR"
    path.write_bytes(raw)
    with WavecarFile(path) as wavecar:
        for indexes in ((-1, 0, 0), (0, 2, 0), (0, 0, 2)):
            with pytest.raises(ValueError, match="out of range"):
                wavecar.coefficients(*indexes)
    with pytest.raises(ValueError, match="closed"):
        wavecar.coefficients(0, 0, 0)

    short_path = tmp_path / "short-coefficients"
    short_path.write_bytes(raw)
    short_wavecar = WavecarFile(short_path)
    os.truncate(short_path, 3 * 128)
    try:
        with pytest.raises(ValueError, match="short"):
            short_wavecar.coefficients(0, 0, 0)
    finally:
        short_wavecar.close()

    with pytest.raises(ValueError, match="decompressed"):
        WavecarFile(tmp_path / "WAVECAR.bz2")
    with pytest.raises(ValueError, match="binary.*compression"):
        _write_wavecar_payload(io.StringIO(), {"format": "vasp-wavecar", "wavecar": Source(data)})


@requires_numpy
def test_wavecar_registration_and_load(tmp_path: pathlib.Path) -> None:
    data = _golden_data()
    path = tmp_path / "WAVECAR"
    path.write_bytes(_golden_bytes(data, 128, numpy.complex64))
    assert httk.core.has_reader_for(str(path))
    assert httk.core.has_writer_for(str(path))

    raw = httk.core.load(str(path), raw=True)
    default = None
    try:
        importlib.import_module("httk.atomistic")
        has_adapter = True
    except ImportError:
        has_adapter = False
    try:
        default = httk.core.load(str(path))
        assert raw["format"] == "vasp-wavecar"
        assert isinstance(raw["wavecar"], WavecarFile)
        if has_adapter:
            assert type(default).__name__ == "PlaneWaveFunctions"
        else:
            assert default["format"] == "vasp-wavecar"
            assert isinstance(default["wavecar"], WavecarFile)
    finally:
        raw["wavecar"].close()
        if default is not None:
            if has_adapter:
                default.close()
            else:
                default["wavecar"].close()


def test_import_and_wavecar_import_without_numpy() -> None:
    src_dir = str(pathlib.Path(__file__).resolve().parent.parent / "src")
    script = (
        "import importlib, sys\n"
        "sys.modules['numpy'] = None\n"
        f"sys.path.insert(0, {src_dir!r})\n"
        "import httk.io\n"
        "assert hasattr(httk.io, 'read_cif')\n"
        "try:\n"
        "    importlib.import_module('httk.io.vasp.wavecar')\n"
        "except ImportError as error:\n"
        "    assert 'httk-io[numpy]' in str(error), error\n"
        "else:\n"
        "    raise AssertionError('wavecar import unexpectedly succeeded')\n"
    )
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
