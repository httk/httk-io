"""Lazy, lexeme-preserving readers for VASP ``XDATCAR`` files."""

import os
import re
from collections import deque
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from types import TracebackType
from typing import Any, Self

from ._text import source_lines

_CONFIGURATION = re.compile(r"^\s*(Direct|Cartesian)\s+configuration\s*=\s*(\d+)", re.IGNORECASE)
_NUMBER = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?$")


@dataclass(frozen=True)
class _Header:
    comment: str
    scale: str
    cell: tuple[tuple[str, str, str], ...]
    symbols: tuple[str, ...] | None
    counts: tuple[int, ...]
    cartesian: bool


class _NumberedLines:
    def __init__(self, lines: Iterator[str]) -> None:
        self._lines = enumerate(lines, 1)
        self._pending: deque[tuple[int, str]] = deque()

    def next(self) -> tuple[int, str]:
        if self._pending:
            return self._pending.popleft()
        return next(self._lines)

    def __iter__(self) -> Self:
        return self

    def __next__(self) -> tuple[int, str]:
        return self.next()

    def take(self, count: int) -> list[tuple[int, str]]:
        return [self.next() for _ in range(count)]

    def push(self, record: tuple[int, str]) -> None:
        self._pending.appendleft(record)


def _tokens(line: str) -> list[str]:
    return line.strip().split()


def _parse_header(records: list[tuple[int, str]]) -> _Header:
    if len(records) not in (6, 7):
        raise ValueError("XDATCAR header is truncated")
    values = [_tokens(line) for _number, line in records]
    if (
        len(values[1]) != 1
        or not _NUMBER.fullmatch(values[1][0])
        or any(len(values[i]) != 3 or any(not _NUMBER.fullmatch(token) for token in values[i]) for i in (2, 3, 4))
        or not values[5]
    ):
        raise ValueError("invalid XDATCAR header")
    if len(records) == 6:
        symbols = None
        counts_tokens = values[5]
    else:
        try:
            tuple(int(token) for token in values[5])
        except ValueError:
            symbols = tuple(values[5])
        else:
            symbols = None
        counts_tokens = values[6]
    try:
        counts = tuple(int(token) for token in counts_tokens)
    except ValueError:
        raise ValueError("invalid XDATCAR atom counts") from None
    if not counts or any(count < 0 for count in counts):
        raise ValueError("invalid XDATCAR atom counts")
    if len(symbols or counts) != len(counts):
        raise ValueError("XDATCAR symbols and counts have different lengths")
    return _Header(
        comment=records[0][1].rstrip("\r\n"),
        scale=values[1][0],
        cell=tuple((values[row][0], values[row][1], values[row][2]) for row in (2, 3, 4)),
        symbols=symbols,
        counts=counts,
        cartesian=False,
    )


def _header_from_source(lines: Iterator[str]) -> _Header:
    return _read_header(_NumberedLines(lines))


def _read_header(reader: _NumberedLines, first: tuple[int, str] | None = None) -> _Header:
    records = ([first] if first is not None else []) + reader.take(6 if first is None else 5)
    try:
        return _parse_header(records)
    except ValueError:
        records.append(reader.next())
        return _parse_header(records)


class XdatcarFile:
    """Re-openable, forward-streaming XDATCAR source.

    Construction checks only that ``filename`` exists.  Header properties scan the
    seven-line header; ``frames`` opens a fresh stream and never caches frames.
    Variable-cell files are identified by a repeated POSCAR-like header and expose
    that header's cell on the following frame.  An incomplete final coordinate block
    is dropped and reported in :attr:`issues` during the full pass. The public
    :attr:`path` property returns the source filename as a string.
    """

    def __init__(self, filename: str | os.PathLike[str]) -> None:
        if not isinstance(filename, (str, os.PathLike)):
            raise TypeError("XDATCAR requires a filesystem filename, not a live stream")
        self._path = Path(filename)
        self._filename = str(filename)
        if not self._path.is_file():
            raise FileNotFoundError(self._path)
        self._closed = False
        self._header: _Header | None = None
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
        """Close the object; scans use short-lived streams and own no handle."""
        self._closed = True

    def __enter__(self) -> Self:
        if self.closed:
            raise ValueError("Cannot enter a closed XDATCAR file.")
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
            raise ValueError("Cannot read a closed XDATCAR file.")

    def _ensure_header(self) -> _Header:
        self._check_open()
        if self._header is None:
            with source_lines(self._path) as (lines, _raw):
                reader = _NumberedLines(iter(lines))
                header = _read_header(reader)
                cartesian = False
                for _number, line in reader:
                    match = _CONFIGURATION.match(line)
                    if match is not None:
                        cartesian = match.group(1).lower() == "cartesian"
                        break
                self._header = _Header(
                    header.comment,
                    header.scale,
                    header.cell,
                    header.symbols,
                    header.counts,
                    cartesian,
                )
        return self._header

    @property
    def comment(self) -> str:
        return self._ensure_header().comment

    @property
    def scale(self) -> str:
        return self._ensure_header().scale

    @property
    def cell(self) -> tuple[tuple[str, str, str], ...]:
        return self._ensure_header().cell

    @property
    def symbols(self) -> tuple[str, ...] | None:
        return self._ensure_header().symbols

    @property
    def counts(self) -> tuple[int, ...]:
        return self._ensure_header().counts

    @property
    def cartesian(self) -> bool:
        return self._ensure_header().cartesian

    @property
    def issues(self) -> tuple[str, ...]:
        self._ensure_full()
        return self._issues or ()

    def _iter_frames(self, *, issues: list[str] | None = None) -> Iterator[Mapping[str, Any]]:
        with source_lines(self._path) as (lines, _raw):
            reader = _NumberedLines(iter(lines))
            header = _read_header(reader)
            variable_cell = False
            current_cell: tuple[tuple[str, str, str], ...] | None = None
            current_scale: str | None = None
            frame_index = 0
            previous_configuration: int | None = None
            atom_count = sum(header.counts)
            while True:
                try:
                    number, line = reader.next()
                except StopIteration:
                    return
                match = _CONFIGURATION.match(line)
                if match is not None:
                    configuration = int(match.group(2))
                    if (
                        previous_configuration is not None
                        and configuration != previous_configuration + 1
                        and issues is not None
                    ):
                        issues.append(
                            f"line {number}: discontinuous configuration index {configuration}; "
                            f"expected {previous_configuration + 1}"
                        )
                    previous_configuration = configuration
                    cartesian = match.group(1).lower() == "cartesian"
                    coordinates: list[tuple[str, str, str]] = []
                    for _ in range(atom_count):
                        try:
                            coordinate_number, coordinate_line = reader.next()
                        except StopIteration:
                            if issues is not None:
                                issues.append(f"line {number}: truncated final coordinate block")
                            return
                        tokens = _tokens(coordinate_line)
                        if len(tokens) < 3:
                            if issues is not None:
                                issues.append(f"line {coordinate_number}: malformed coordinate row")
                            return
                        coordinates.append((tokens[0], tokens[1], tokens[2]))
                    frame = {
                        "index": frame_index,
                        "cell": current_cell if variable_cell else None,
                        "coords": tuple(coordinates),
                        "cartesian": cartesian,
                        "scale": current_scale if variable_cell else None,
                    }
                    frame_index += 1
                    yield frame
                    continue
                if not line.strip():
                    continue
                # NPT XDATCAR repeats the complete seven-line header before a marker.
                try:
                    repeated_header = _read_header(reader, (number, line))
                except (StopIteration, ValueError):
                    if issues is not None:
                        issues.append(f"line {number}: malformed XDATCAR header or unexpected text")
                    return
                if repeated_header.symbols != header.symbols or repeated_header.counts != header.counts:
                    if issues is not None:
                        issues.append(f"line {number}: repeated XDATCAR symbols/counts differ from initial header")
                    return
                try:
                    marker_number, marker_line = reader.next()
                    while not marker_line.strip():
                        marker_number, marker_line = reader.next()
                except StopIteration:
                    if issues is not None:
                        issues.append(f"line {number}: repeated header has no following configuration marker")
                    return
                if _CONFIGURATION.match(marker_line) is None:
                    if issues is not None:
                        issues.append(f"line {marker_number}: repeated header not followed by configuration marker")
                    return
                reader.push((marker_number, marker_line))
                variable_cell = True
                current_cell = repeated_header.cell
                current_scale = repeated_header.scale
                atom_count = sum(repeated_header.counts)

    def frames(self) -> Iterator[Mapping[str, Any]]:
        """Yield complete frames in file order without retaining them."""
        self._check_open()
        issues: list[str] = []
        yield from self._iter_frames(issues=issues)
        self._issues = tuple(issues)

    def _ensure_full(self) -> None:
        self._check_open()
        if self._nframes is not None:
            return
        issues: list[str] = []
        count = 0
        for _frame in self._iter_frames(issues=issues):
            count += 1
        self._nframes = count
        self._issues = tuple(issues)

    @property
    def nframes(self) -> int:
        self._ensure_full()
        return self._nframes or 0


def read_xdatcar(source: str | os.PathLike[str]) -> dict[str, Any]:
    """Read an XDATCAR path into a lazy neutral ``vasp-xdatcar`` payload."""
    return {"format": "vasp-xdatcar", "xdatcar": XdatcarFile(source)}
