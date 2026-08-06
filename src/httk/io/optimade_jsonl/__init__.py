"""Streaming OPTIMADE partial-data JSON Lines trajectories.

The stable ``httk-trajectory-jsonl`` 0.1 format is one JSON object per line.
The first line is the header::

    {
      "optimade-partial-data": {"format": "1.2.0"},
      "layout": "dense",
      "x-httk-trajectory": {
        "format": "httk-trajectory-jsonl", "version": "0.1",
        "species": [{"name": "Si", "chemical_symbols": ["Si"],
                     "concentration": [1.0]}],
        "species_at_sites": ["Si"],
        "constant_cell": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0]],
        "nframes": 3,
        "observable_names": ["energy"],
        "reference_frames": [0],
        "line_schema": {"...": "see below"}
      }
    }

Each following line is a frame object with exactly these semantic members::

    {"index": 0,
     "fractional_site_positions": [[0.0, 0.0, 0.0]],
     "observables": {"energy": -1.25}}

``index`` is the zero-based frame number.  ``fractional_site_positions`` is an
N-by-3 array of float64 presentation values.  ``observables`` contains every
name declared by the header, with ``null`` allowed.  If ``constant_cell`` is
null, each frame additionally contains ``lattice_vectors`` as a 3-by-3 array
of float64 presentation values.  If it is non-null, a frame may include that
member only when it is equal to the declared constant cell.  The writer emits
the compact form and the reader accepts either form.  No exact-token channel
exists: JSON numbers intentionally follow the float64 presentation, like the
numeric layer.

This is a holding/container format, not a database representation.  It follows
the OPTIMADE partial-data JSON Lines framing (OPTIMADE 1.2, dense layout), but
uses the ``x-httk-trajectory`` member because one line carries a complete frame
and several properties.  A binary framing variant is deliberately deferred
pending a separate design decision.
"""

import bz2
import gzip
import json
import lzma
import math
import os
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from types import TracebackType
from typing import Any, Self, TextIO

from httk.core import TextstreamFileView
from httk.core.datastream.compression import split_compression_suffix

FORMAT = "httk-trajectory-jsonl"
VERSION = "0.1"
_PARTIAL_FORMAT = "1.2.0"


def _line_schema(constant_cell: Any) -> dict[str, Any]:
    required = ["index", "fractional_site_positions", "observables"]
    properties: dict[str, Any] = {
        "index": {"type": "integer", "minimum": 0},
        "fractional_site_positions": {"type": "array", "items": {"type": "array"}},
        "observables": {"type": "object"},
        "lattice_vectors": {"type": "array", "items": {"type": "array"}},
    }
    if constant_cell is None:
        required.append("lattice_vectors")
    return {"type": "object", "required": required, "properties": properties}


def _as_float(value: Any, *, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field} must contain JSON numbers")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{field} must contain finite JSON numbers")
    return result


def _matrix(value: Any, rows: int, columns: int, *, field: str) -> list[list[float]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != rows:
        raise ValueError(f"{field} must be a {rows}x{columns} array")
    result = []
    for row in value:
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or len(row) != columns:
            raise ValueError(f"{field} must be a {rows}x{columns} array")
        result.append([_as_float(item, field=field) for item in row])
    return result


def _same_matrix(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> bool:
    return all(a == b for row_a, row_b in zip(left, right) for a, b in zip(row_a, row_b))


def _header_info(header: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    if not isinstance(header, Mapping):
        raise TypeError("trajectory JSONL header must be a mapping")
    if "x-httk-trajectory" in header:
        partial = dict(header)
        info = partial.get("x-httk-trajectory")
        if not isinstance(info, Mapping):
            raise ValueError("x-httk-trajectory must be a mapping")
        info = dict(info)
    else:
        info = dict(header)
        info.setdefault("format", FORMAT)
        info.setdefault("version", VERSION)
        partial = {
            "optimade-partial-data": {"format": _PARTIAL_FORMAT},
            "layout": "dense",
            "x-httk-trajectory": info,
        }
    partial_data = partial.get("optimade-partial-data")
    if partial_data != {"format": _PARTIAL_FORMAT} and (
        not isinstance(partial_data, Mapping) or partial_data.get("format") != _PARTIAL_FORMAT
    ):
        raise ValueError(f"optimade-partial-data.format must be {_PARTIAL_FORMAT!r}")
    if partial.get("layout") != "dense":
        raise ValueError("trajectory JSONL requires the OPTIMADE dense layout")
    if info.get("format") != FORMAT or info.get("version") != VERSION:
        raise ValueError(f"x-httk-trajectory must declare format={FORMAT!r}, version={VERSION!r}")
    for key in ("species", "species_at_sites", "observable_names", "reference_frames"):
        if key not in info:
            raise ValueError(f"trajectory JSONL header is missing {key!r}")
    if not isinstance(info["species"], Sequence) or isinstance(info["species"], (str, bytes)):
        raise ValueError("trajectory JSONL header species must be an array")
    if not isinstance(info["species_at_sites"], Sequence) or isinstance(info["species_at_sites"], (str, bytes)):
        raise ValueError("trajectory JSONL header species_at_sites must be an array")
    if len(info["species_at_sites"]) < 1 or any(not isinstance(name, str) for name in info["species_at_sites"]):
        raise ValueError("trajectory JSONL header must describe at least one site")
    if not isinstance(info["observable_names"], Sequence) or isinstance(info["observable_names"], (str, bytes)):
        raise ValueError("trajectory JSONL header observable_names must be an array")
    if any(not isinstance(name, str) for name in info["observable_names"]):
        raise ValueError("trajectory JSONL observable_names must contain strings")
    if len(set(info["observable_names"])) != len(info["observable_names"]):
        raise ValueError("trajectory JSONL observable_names must be unique")
    cell = info.get("constant_cell")
    if cell is not None:
        info["constant_cell"] = _matrix(cell, 3, 3, field="constant_cell")
    nframes = info.get("nframes")
    if nframes is not None and (not isinstance(nframes, int) or isinstance(nframes, bool) or nframes < 1):
        raise ValueError("trajectory JSONL nframes must be a positive integer or null")
    references = info["reference_frames"]
    if references is not None and (
        not isinstance(references, Sequence)
        or isinstance(references, (str, bytes))
        or any(not isinstance(value, int) or isinstance(value, bool) or value < 0 for value in references)
    ):
        raise ValueError("trajectory JSONL reference_frames must be an array of non-negative integers or null")
    if nframes is not None and references is not None and any(value >= nframes for value in references):
        raise ValueError("trajectory JSONL reference_frames must be within nframes")
    info["line_schema"] = info.get("line_schema", _line_schema(info.get("constant_cell")))
    partial["x-httk-trajectory"] = info
    return partial, info


def _frame(frame: Mapping[str, Any], info: Mapping[str, Any], expected_index: int) -> dict[str, Any]:
    if not isinstance(frame, Mapping):
        raise TypeError("trajectory JSONL frames must be mappings")
    index = frame.get("index", expected_index)
    if index != expected_index or isinstance(index, bool) or not isinstance(index, int) or index < 0:
        raise ValueError(f"trajectory JSONL frame index must be {expected_index}")
    positions = _matrix(
        frame.get("fractional_site_positions"),
        len(info["species_at_sites"]),
        3,
        field="fractional_site_positions",
    )
    names = tuple(info["observable_names"])
    observables = frame.get("observables", {})
    if not isinstance(observables, Mapping) or set(observables) != set(names):
        raise ValueError("trajectory JSONL frame observables must match header observable_names")
    result: dict[str, Any] = {
        "index": index,
        "fractional_site_positions": positions,
        "observables": {name: _json_numbers(observables[name], field=f"observable {name!r}") for name in names},
    }
    constant = info.get("constant_cell")
    cell = frame.get("lattice_vectors")
    if constant is None:
        if cell is None:
            raise ValueError("variable-cell trajectory frame is missing lattice_vectors")
        result["lattice_vectors"] = _matrix(cell, 3, 3, field="lattice_vectors")
    elif cell is not None and not _same_matrix(constant, _matrix(cell, 3, 3, field="lattice_vectors")):
        raise ValueError("trajectory frame lattice_vectors differs from header constant_cell")
    return result


def _json_numbers(value: Any, *, field: str) -> Any:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return value
    if isinstance(value, (int, float)):
        return _as_float(value, field=field)
    if isinstance(value, Mapping):
        return {str(key): _json_numbers(item, field=field) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [_json_numbers(item, field=field) for item in value]
    raise ValueError(f"{field} contains a value that is not JSON-serializable")


class TrajectoryJsonlFile:
    """A closed-by-default, re-openable, lazy trajectory JSONL reader.

    The public :attr:`path` property returns the source filename as a string.
    """

    def __init__(self, filename: str | os.PathLike[str]) -> None:
        if not isinstance(filename, (str, os.PathLike)):
            raise TypeError("trajectory JSONL requires a filesystem filename")
        self._path = Path(filename)
        self._filename = str(filename)
        if not self._path.is_file():
            raise FileNotFoundError(self._path)
        self._closed = False
        self._header: dict[str, Any] | None = None
        self._nframes: int | None = None
        self._issues: tuple[str, ...] | None = None

    @property
    def path(self) -> str:
        """Return the source filename used to construct this lazy reader."""
        self._check_open()
        return self._filename

    @property
    def closed(self) -> bool:
        return self._closed

    def close(self) -> None:
        self._closed = True

    def __enter__(self) -> Self:
        self._check_open()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def _check_open(self) -> None:
        if self.closed:
            raise ValueError("Cannot read a closed trajectory JSONL file.")

    def _read_header(self) -> tuple[dict[str, Any], dict[str, Any]]:
        self._check_open()
        if self._header is None:
            with TextstreamFileView(self._path) as stream:
                line = stream.readline()
            if not line:
                raise ValueError("trajectory JSONL file is empty")
            try:
                raw = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"invalid trajectory JSONL header: {exc.msg}") from exc
            self._header, _info = _header_info(raw)
        return self._header, self._header["x-httk-trajectory"]

    @property
    def header(self) -> Mapping[str, Any]:
        return self._read_header()[0]

    @property
    def issues(self) -> tuple[str, ...]:
        self._ensure_full()
        return self._issues or ()

    def _iter_frames(self, issues: list[str] | None = None) -> Iterator[Mapping[str, Any]]:
        _header, info = self._read_header()
        with TextstreamFileView(self._path) as stream:
            stream.readline()
            index = 0
            for line_number, line in enumerate(stream, 2):
                if not line.strip():
                    continue
                try:
                    raw = json.loads(line)
                    if not isinstance(raw, Mapping):
                        raise ValueError("frame line must be a JSON object")
                    yield _frame(raw, info, index)
                except (json.JSONDecodeError, TypeError, ValueError) as exc:
                    if issues is not None:
                        issues.append(f"line {line_number}: {exc}")
                        return
                    raise ValueError(f"line {line_number}: invalid trajectory frame: {exc}") from exc
                index += 1

    def frames(self) -> Iterator[Mapping[str, Any]]:
        self._check_open()
        issues: list[str] = []
        count = 0
        for frame in self._iter_frames(issues):
            count += 1
            yield frame
        declared = self._read_header()[1].get("nframes")
        if declared is not None and declared != count:
            issues.append(f"header nframes={declared} does not match {count} frame lines")
        self._nframes = count if declared is None else declared
        self._issues = tuple(issues)

    def _ensure_full(self) -> None:
        self._check_open()
        if self._nframes is not None and self._issues is not None:
            return
        issues: list[str] = []
        count = sum(1 for _frame in self._iter_frames(issues))
        declared = self._read_header()[1].get("nframes")
        if declared is not None and declared != count:
            issues.append(f"header nframes={declared} does not match {count} frame lines")
        self._nframes = count if declared is None else declared
        self._issues = tuple(issues)

    @property
    def nframes(self) -> int:
        self._check_open()
        if self._nframes is None:
            declared = self._read_header()[1].get("nframes")
            if declared is not None:
                self._nframes = declared
            else:
                self._ensure_full()
        return self._nframes or 0

    def frame(self, i: int) -> Mapping[str, Any]:
        self._check_open()
        index = i if i >= 0 else self.nframes + i
        if index < 0:
            raise IndexError(f"trajectory frame index out of range: {i}")
        try:
            return next(frame for frame in self._iter_frames() if frame["index"] == index)
        except StopIteration:
            raise IndexError(f"trajectory frame index out of range: {i}") from None


def write_trajectory_jsonl(
    destination: str | os.PathLike[str] | TextIO,
    header: Mapping[str, Any],
    frames: Iterable[Mapping[str, Any]],
) -> None:
    """Write a validated trajectory header and frame iterator without buffering frames."""
    full_header, info = _header_info(header)
    close = False
    stream: TextIO
    if isinstance(destination, (str, os.PathLike)):
        name = os.fspath(destination)
        _inner, codec = split_compression_suffix(os.path.basename(name))
        if codec is None:
            stream = open(name, "w", encoding="utf-8", newline="\n")  # noqa: SIM115
        elif codec.name == "gzip":
            stream = gzip.open(name, "wt", encoding="utf-8", newline="\n")  # noqa: SIM115
        elif codec.name == "bzip2":
            stream = bz2.open(name, "wt", encoding="utf-8", newline="\n")  # noqa: SIM115
        elif codec.name == "xz":
            stream = lzma.open(name, "wt", encoding="utf-8", newline="\n")  # noqa: SIM115
        elif codec.name == "lzma":
            stream = lzma.open(  # noqa: SIM115
                name, "wt", encoding="utf-8", format=lzma.FORMAT_ALONE, newline="\n"
            )
        else:
            raise ValueError(f"trajectory JSONL cannot write compression codec {codec.name!r}")
        close = True
    else:
        stream = destination
    try:
        stream.write(json.dumps(full_header, separators=(",", ":"), allow_nan=False) + "\n")
        count = 0
        for count, frame in enumerate(frames, 1):
            stream.write(json.dumps(_frame(frame, info, count - 1), separators=(",", ":"), allow_nan=False) + "\n")
        declared = info.get("nframes")
        if declared is not None and declared != count:
            raise ValueError(f"trajectory JSONL header nframes={declared} does not match {count} frames")
    finally:
        if close:
            stream.close()


def read_trajectory_jsonl(source: str | os.PathLike[str]) -> dict[str, Any]:
    """Read a trajectory JSONL path into a lazy neutral payload."""
    return {"format": FORMAT, "trajectory_jsonl": TrajectoryJsonlFile(source)}


def _write_trajectory_jsonl_payload(destination: TextIO, payload: Mapping[str, Any]) -> None:
    if payload.get("format") != FORMAT:
        raise ValueError(f"expected {FORMAT!r} payload")
    write_trajectory_jsonl(destination, payload["header"], payload["frames"])


__all__ = ["FORMAT", "VERSION", "TrajectoryJsonlFile", "read_trajectory_jsonl", "write_trajectory_jsonl"]
