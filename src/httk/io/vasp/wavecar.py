"""Neutral, numpy-backed VASP WAVECAR binary I/O.

The ``wavecar`` value returned by :func:`read_wavecar` is a WavecarSource:
an object exposing ``nspins``, ``nkpts``, ``nbands``, ``encut``, ``cell``,
``kpoints``, ``eigenvalues``, ``occupations``, ``nplanewaves``,
``double_precision``, and ``coefficients(spin, kpt, band)``. ``WavecarFile``
also exposes its ``record_length``. The gamma-half orientation is not stored
in a WAVECAR; ``gamma_half`` is a pass-through hint because the consumer's
default (``"x"``) cannot be inferred from the data.
The atomistic phase may provide the same contract for in-memory data without
depending on this module's concrete class.
"""

import io
import os
from collections.abc import Mapping
from types import TracebackType
from typing import Any, BinaryIO, Self

from httk.core.datastream.compression import split_compression_suffix

try:
    import numpy
except ImportError:
    raise ImportError("httk.io.vasp.wavecar requires numpy; install httk-io[numpy]") from None


_HEADER_FLOATS = 12
_HEADER_BYTES = _HEADER_FLOATS * numpy.dtype(numpy.float64).itemsize


def _read_floats(file: BinaryIO, count: int, description: str) -> Any:
    values = numpy.fromfile(file, dtype=numpy.float64, count=count)
    if len(values) != count:
        raise ValueError(f"WAVECAR is truncated while reading {description}.")
    return values


def _positive_integer(value: Any, name: str) -> int:
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"WAVECAR {name} must be a positive integer, got {value!r}.") from None
    if not numpy.isfinite(number) or not number.is_integer() or number <= 0:
        raise ValueError(f"WAVECAR {name} must be a positive integer, got {value!r}.")
    return int(number)


def _index(value: int, size: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, numpy.integer)) or not 0 <= value < size:
        raise ValueError(f"WAVECAR {name} index {value!r} is out of range [0, {size}).")
    return int(value)


class WavecarFile:
    """Open a VASP WAVECAR and eagerly read its small metadata headers.

    Public indices are all zero-based. ``cell`` has shape ``(3, 3)``,
    ``kpoints`` has shape ``(nkpts, 3)``, the eigenvalue and occupation arrays
    have shape ``(nspins, nkpts, nbands)``, and ``nplanewaves`` has shape
    ``(nkpts,)``. Coefficients are read afresh for each call.
    """

    nspins: int
    nkpts: int
    nbands: int
    encut: float
    cell: Any
    kpoints: Any
    eigenvalues: Any
    occupations: Any
    nplanewaves: Any
    double_precision: bool
    record_length: int

    def __init__(self, filename: str | os.PathLike[str], *, double_precision: bool | None = None) -> None:
        if double_precision is not None and not isinstance(double_precision, bool):
            raise ValueError("WAVECAR double_precision must be True, False, or None.")
        path = os.fspath(filename)
        basename = os.fsdecode(os.path.basename(path))
        _inner, compression = split_compression_suffix(basename)
        if compression is not None:
            raise ValueError(
                "WAVECAR files must be decompressed on disk first; random access cannot stream decompression."
            )

        file = open(path, "rb")  # noqa: SIM115 - the handle is owned until close() or context exit.
        self._file: BinaryIO | None = file
        try:
            self._read_metadata(double_precision)
        except Exception:
            file.close()
            self._file = None
            raise

    def _read_metadata(self, precision_override: bool | None) -> None:
        file = self._file
        if file is None:
            raise ValueError("Cannot read metadata from a closed WAVECAR file.")
        record_zero = _read_floats(file, 3, "record 0")
        self.record_length = _positive_integer(record_zero[0], "record length")
        # Record 1 has 12 float64 values and defines the minimum record length.
        if self.record_length < _HEADER_BYTES:
            raise ValueError(
                f"WAVECAR record length {self.record_length} is smaller than the {_HEADER_BYTES}-byte record 1 header."
            )

        rtag = record_zero[2]
        if precision_override is None:
            if rtag == 45200:
                self.double_precision = False
            elif rtag == 45210:
                self.double_precision = True
            else:
                raise ValueError(f"Unknown WAVECAR RTAG value {rtag!r}.")
        else:
            self.double_precision = precision_override
        self._coefficient_dtype = numpy.complex128 if self.double_precision else numpy.complex64

        file.seek(self.record_length)
        record_one = _read_floats(file, _HEADER_FLOATS, "record 1")
        self.nkpts = _positive_integer(record_one[0], "number of k-points")
        self.nbands = _positive_integer(record_one[1], "number of bands")
        self.encut = float(record_one[2])
        self.cell = numpy.array(record_one[3:12], dtype=numpy.float64).reshape(3, 3)
        self.nspins = _positive_integer(record_zero[1], "number of spins")

        header_bytes = (4 + 3 * self.nbands) * numpy.dtype(numpy.float64).itemsize
        if self.record_length < header_bytes:
            raise ValueError(
                f"WAVECAR record length {self.record_length} is too small for {self.nbands} band header entries."
            )
        record_count = 2 + self.nspins * self.nkpts * (self.nbands + 1)
        file_size = os.fstat(file.fileno()).st_size
        expected_size = record_count * self.record_length
        if file_size < expected_size:
            raise ValueError(f"WAVECAR is truncated: expected at least {expected_size} bytes, found {file_size}.")

        self.kpoints = numpy.zeros((self.nkpts, 3), dtype=numpy.float64)
        self.eigenvalues = numpy.zeros((self.nspins, self.nkpts, self.nbands), dtype=numpy.float64)
        self.occupations = numpy.zeros_like(self.eigenvalues)
        self.nplanewaves = numpy.zeros(self.nkpts, dtype=numpy.int64)

        for spin in range(self.nspins):
            for kpt in range(self.nkpts):
                header_record = 2 + spin * self.nkpts * (self.nbands + 1) + kpt * (self.nbands + 1)
                file.seek(header_record * self.record_length)
                header = _read_floats(file, 4 + 3 * self.nbands, f"spin {spin}, k-point {kpt} header")
                nplanewaves = _positive_integer(header[0], f"spin {spin}, k-point {kpt} plane-wave count")
                kpoint = header[1:4]
                if spin == 0:
                    self.nplanewaves[kpt] = nplanewaves
                    self.kpoints[kpt] = kpoint
                elif nplanewaves != self.nplanewaves[kpt] or not numpy.array_equal(kpoint, self.kpoints[kpt]):
                    raise ValueError(f"WAVECAR k-point {kpt} metadata differs between spin channels.")
                self.eigenvalues[spin, kpt] = header[4::3]
                self.occupations[spin, kpt] = header[6::3]

    @property
    def closed(self) -> bool:
        """Whether the underlying WAVECAR file is closed."""
        return self._file is None or self._file.closed

    def close(self) -> None:
        """Close the WAVECAR file."""
        if self._file is not None:
            self._file.close()

    def __enter__(self) -> Self:
        if self.closed:
            raise ValueError("Cannot enter a closed WAVECAR file.")
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def coefficients(self, spin: int, kpt: int, band: int) -> Any:
        """Read one zero-based spin/k-point/band coefficient vector."""
        if self.closed:
            raise ValueError("Cannot read coefficients from a closed WAVECAR file.")
        spin = _index(spin, self.nspins, "spin")
        kpt = _index(kpt, self.nkpts, "k-point")
        band = _index(band, self.nbands, "band")
        header_record = 2 + spin * self.nkpts * (self.nbands + 1) + kpt * (self.nbands + 1)
        coefficient_record = header_record + 1 + band
        file = self._file
        if file is None:
            raise ValueError("Cannot read coefficients from a closed WAVECAR file.")
        file.seek(coefficient_record * self.record_length)
        values = numpy.fromfile(file, dtype=self._coefficient_dtype, count=int(self.nplanewaves[kpt]))
        if len(values) != self.nplanewaves[kpt]:
            raise ValueError(
                f"WAVECAR coefficient record is short: expected {self.nplanewaves[kpt]} values, found {len(values)}."
            )
        return values


def read_wavecar(
    source: str | os.PathLike[str], *, double_precision: bool | None = None, gamma_half: str | None = None
) -> dict[str, Any]:
    """Read a WAVECAR path into a neutral ``vasp-wavecar`` payload."""
    if gamma_half not in (None, "x", "z"):
        raise ValueError("gamma_half must be None, 'x', or 'z'.")
    return {
        "format": "vasp-wavecar",
        "wavecar": WavecarFile(source, double_precision=double_precision),
        "gamma_half": gamma_half,
    }


def _write_wavecar_payload(destination: str | os.PathLike[str], data: Mapping[str, Any], **kwargs: Any) -> None:
    """Write a neutral WAVECAR payload through the core writer registry."""
    if kwargs:
        unexpected = ", ".join(sorted(kwargs))
        raise TypeError(f"Unexpected WAVECAR writer keyword argument(s): {unexpected}.")
    if isinstance(destination, io.TextIOBase):
        raise ValueError("WAVECAR is binary and cannot be written through compression.")
    if not isinstance(destination, (str, os.PathLike)):
        raise ValueError("WAVECAR destination must be a filesystem path; WAVECAR is binary.")
    if not isinstance(data, Mapping):
        raise TypeError("WAVECAR data must be a neutral payload mapping.")
    if data.get("format", "vasp-wavecar") != "vasp-wavecar":
        raise ValueError("WAVECAR payload format must be 'vasp-wavecar'.")
    source = data["wavecar"]
    nspins = _positive_integer(source.nspins, "number of spins")
    nkpts = _positive_integer(source.nkpts, "number of k-points")
    nbands = _positive_integer(source.nbands, "number of bands")
    if not isinstance(source.double_precision, bool):
        raise ValueError("WAVECAR source double_precision must be a boolean.")
    nplanewaves = numpy.asarray(source.nplanewaves)
    if nplanewaves.shape != (nkpts,):
        raise ValueError(f"WAVECAR source nplanewaves must have shape ({nkpts},).")
    if not numpy.all(numpy.isfinite(nplanewaves)) or not numpy.all(nplanewaves == nplanewaves.astype(numpy.int64)):
        raise ValueError("WAVECAR source nplanewaves must contain integers.")
    nplanewaves = nplanewaves.astype(numpy.int64)
    if numpy.any(nplanewaves <= 0):
        raise ValueError("WAVECAR source nplanewaves must be positive.")
    kpoints = numpy.asarray(source.kpoints, dtype=numpy.float64)
    cell = numpy.asarray(source.cell, dtype=numpy.float64)
    eigenvalues = numpy.asarray(source.eigenvalues, dtype=numpy.float64)
    occupations = numpy.asarray(source.occupations, dtype=numpy.float64)
    if kpoints.shape != (nkpts, 3) or cell.shape != (3, 3):
        raise ValueError("WAVECAR source kpoints and cell must have shapes (nkpts, 3) and (3, 3).")
    expected_shape = (nspins, nkpts, nbands)
    if eigenvalues.shape != expected_shape or occupations.shape != expected_shape:
        raise ValueError(f"WAVECAR source eigenvalues and occupations must have shape {expected_shape}.")

    coefficient_dtype = numpy.complex128 if source.double_precision else numpy.complex64
    coefficient_itemsize = numpy.dtype(coefficient_dtype).itemsize
    record_length = max(
        int(numpy.max(nplanewaves)) * coefficient_itemsize,
        (4 + 3 * nbands) * 8,
        12 * 8,
    )
    record_length = ((record_length + coefficient_itemsize - 1) // coefficient_itemsize) * coefficient_itemsize
    record = numpy.zeros(record_length // 8, dtype=numpy.float64)

    with open(destination, "wb") as file:
        record[:3] = (record_length, nspins, 45210 if source.double_precision else 45200)
        record.tofile(file)
        record.fill(0)
        record[:12] = numpy.concatenate(([nkpts, nbands, float(source.encut)], cell.reshape(-1)))
        record.tofile(file)
        for spin in range(nspins):
            for kpt in range(nkpts):
                record.fill(0)
                record[0] = nplanewaves[kpt]
                record[1:4] = kpoints[kpt]
                for band in range(nbands):
                    record[4 + 3 * band] = eigenvalues[spin, kpt, band]
                    record[6 + 3 * band] = occupations[spin, kpt, band]
                record.tofile(file)
                for band in range(nbands):
                    coefficients = numpy.asarray(source.coefficients(spin, kpt, band))
                    if coefficients.ndim != 1 or len(coefficients) != nplanewaves[kpt]:
                        raise ValueError(
                            f"WAVECAR source coefficients({spin}, {kpt}, {band}) must contain "
                            f"exactly {nplanewaves[kpt]} values."
                        )
                    coefficients = coefficients.astype(coefficient_dtype, copy=False)
                    file.write(coefficients.tobytes(order="C"))
                    file.write(b"\0" * (record_length - len(coefficients) * coefficient_itemsize))


def write_wavecar(destination: str | os.PathLike[str], payload: Mapping[str, Any]) -> None:
    """Write a neutral ``vasp-wavecar`` payload to a binary path."""
    _write_wavecar_payload(destination, payload)
