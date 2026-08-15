"""Lazy, streaming access to the bounded metadata and final results of OUTCAR."""

import os
import re
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from itertools import islice
from types import MappingProxyType, TracebackType
from typing import Any, Self

from ._text import source_lines

# Ported from old/httk/src/httk/iface/vasp_if.py:get_computation_info.
vasp_xc_tags = {
    "91": "Perdew - Wang 91 (PW-91)",
    "PE": "Perdew-Burke-Ernzerhof (PBE)",
    "AM": "AM05",
    "HL": "Hedin-Lundqvist",
    "CA": "Ceperley-Alder",
    "PZ": "Ceperley-Alder, parametrization of Perdew-Zunger",
    "WI": "Wigner",
    "RP": "revised Perdew-Burke-Ernzerhof (RPBE) with Pade Approximation",
    "RE": "revPBE",
    "VW": "Vosko-Wilk-Nusair (VWN)",
    "B3": "B3LYP, where LDA part is with VWN3-correlation",
    "B5": "B3LYP, where LDA part is with VWN5-correlation",
    "BF": "BEEF, xc (with libbeef)",
    "CO": "no exchange-correlation",
    "PS": "Perdew-Burke-Ernzerhof revised for solids (PBEsol)",
    "LIBXC": "LDA or GGA from Libxc",
    "LI": "LDA or GGA from Libxc",
    "OR": "optPBE",
    "BO": "optB88",
    "MK": "optB86b",
    "RA": "new RPA Perdew Wang",
    "03": "range-separated ACFDT (LDA - sr RPA) mu = 0.3Å",
    "05": "range-separated ACFDT (LDA - sr RPA) mu = 0.5Å",
    "10": "range-separated ACFDT (LDA - sr RPA) mu = 1.0Å",
    "20": "range-separated ACFDT (LDA - sr RPA) mu = 2.0Å",
    "PL": "new RPA+ Perdew Wang",
}

_PARAMETERS = (
    "ENCUT",
    "NKPTS",
    "POTIM",
    "NSW",
    "IBRION",
    "ISIF",
    "NELM",
    "EDIFF",
    "EDIFFG",
    "ISMEAR",
    "SIGMA",
    "ISPIN",
    "GGA",
    "LEXCH",
)
_NUMERIC_PARAMETERS = frozenset(_PARAMETERS[:-2])
_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
_IONIC_MARKER = re.compile(r"\bIteration\s+\d+\s*\(")
_FRAME_MARKER = re.compile(r"\bIteration\s+\d+\s*\(\s*1\s*\)")
_VERSION_TOKEN = re.compile(r"^vasp\.")
_TITEL = re.compile(r"^\s*TITEL\s*=\s*(.*?)\s*$")
_COMPLETION = "General timing and accounting informations"
_ELASTIC_HEADINGS = (
    "TOTAL ELASTIC MODULI",
    "SYMMETRIZED ELASTIC MODULI",
    "ELASTIC MODULI CONTR FROM IONIC RELAXATION",
)


@dataclass(frozen=True)
class FinalEnergies:
    """Store energy lexemes from the last complete or partial energy block.

    :param free_energy: Free-energy lexeme from the block, when present.
    :param energy_without_entropy: Energy-without-entropy lexeme, when present.
    :param energy_sigma0: Sigma-zero energy lexeme, when present.
    :param final: Whether the block was complete when parsing ended.
    """

    free_energy: str | None
    energy_without_entropy: str | None
    energy_sigma0: str | None
    final: bool


@dataclass(frozen=True)
class OutcarFrame:
    """Store one complete ionic-step snapshot with VASP numeric lexemes unchanged.

    :param index: Zero-based ionic-step index.
    :param cell: Lattice-vector lexemes, when the step contains a cell.
    :param positions: Position lexemes, when the step contains positions.
    :param forces: Force lexemes, when the step contains forces.
    :param stress_kbar: Six stress lexemes in VASP order, when present.
    :param free_energy: Free-energy lexeme, when present.
    :param energy_without_entropy: Energy-without-entropy lexeme, when present.
    :param energy_sigma0: Sigma-zero energy lexeme, when present.
    :param temperature: Temperature lexeme, when present.
    """

    index: int
    cell: tuple[tuple[str, str, str], ...] | None
    positions: tuple[tuple[str, str, str], ...] | None
    forces: tuple[tuple[str, str, str], ...] | None
    stress_kbar: tuple[str, ...] | None
    free_energy: str | None
    energy_without_entropy: str | None
    energy_sigma0: str | None
    temperature: str | None

    @staticmethod
    def _floats(rows: tuple[tuple[str, str, str], ...] | None) -> tuple[tuple[float, float, float], ...] | None:
        return None if rows is None else tuple((float(row[0]), float(row[1]), float(row[2])) for row in rows)

    def cell_floats(self) -> tuple[tuple[float, float, float], ...] | None:
        """Convert the cell lexemes to floating-point values when present.

        :return: Converted cell values, or ``None`` when the frame has no cell.
        """
        return self._floats(self.cell)

    def positions_floats(self) -> tuple[tuple[float, float, float], ...] | None:
        """Convert the position lexemes to floating-point values when present.

        :return: Converted position values, or ``None`` when the frame has no positions.
        """
        return self._floats(self.positions)

    def forces_floats(self) -> tuple[tuple[float, float, float], ...] | None:
        """Convert the force lexemes to floating-point values when present.

        :return: Converted force values, or ``None`` when the frame has no forces.
        """
        return self._floats(self.forces)

    def stress_gpa_voigt(self) -> tuple[float, ...] | None:
        """Convert stress to tensile-positive GPa Voigt order ``xx, yy, zz, yz, xz, xy``.

        The conversion multiplies kbar by ``0.1``, reverses VASP's
        compressive-positive sign, and reorders the shear components.

        :return: Stress in tensile-positive GPa Voigt order, or ``None`` when absent.
        """
        if self.stress_kbar is None:
            return None
        values = tuple(float(token) for token in self.stress_kbar)
        # VASP writes [xx yy zz xy yz zx] in kB, with compressive-positive signs.
        return tuple(-values[index] * 0.1 for index in (0, 1, 2, 4, 5, 3))


@dataclass(frozen=True)
class ElasticModuliBlock:
    """Store one six-by-six elastic-moduli table.

    :param heading: Heading identifying the table in the source.
    :param rows: Table rows with their source numeric lexemes.
    """

    heading: str
    rows: tuple[tuple[str, ...], ...]


@dataclass(frozen=True)
class _Prologue:
    version_string: str
    version_numbers: tuple[int, ...]
    nions: int | None
    ions_per_type: tuple[int, ...] | None
    parameters: Mapping[str, str]
    xc: str | None
    potcar_titles: tuple[str, ...]
    issues: tuple[str, ...]
    parameter_lines: Mapping[str, int]


@dataclass(frozen=True)
class _FullPass:
    final_energies: FinalEnergies
    completed: bool
    completion_evidence: tuple[str, ...]
    issues: tuple[str, ...]
    nframes: int
    last_frame: OutcarFrame | None
    stresses: tuple[tuple[str, ...], ...]
    elastic_moduli: tuple[ElasticModuliBlock, ...]
    magnetization: tuple[float, ...] | None
    noncollinear_magnetization: bool


def _issue(issues: list[str], line: int, message: str) -> None:
    issues.append(f"line {line}: {message}")


def _version_numbers(version: str) -> tuple[int, ...]:
    if not version:
        return ()
    numbers: list[int] = []
    for part in version.removeprefix("vasp.").split("."):
        if not part.isdigit():
            break
        numbers.append(int(part))
    return tuple(numbers)


def _number(token: str) -> bool:
    return re.fullmatch(_NUMBER, token) is not None


def _stress_row(line: str) -> tuple[str, ...] | None:
    match = re.search(r"\bin kB\b(.*)$", line)
    if match is None:
        return None
    tokens = match.group(1).split()
    return tuple(tokens[:6]) if len(tokens) >= 6 else None


class _EnergyTracker:
    def __init__(self) -> None:
        self.seen = False
        self.active = False
        self.incomplete = False
        self.block_line = 0
        self.truncated_lines: list[int] = []
        self.free_energy: str | None = None
        self.energy_without_entropy: str | None = None
        self.energy_sigma0: str | None = None

    def feed(self, line: str, lineno: int) -> None:
        if "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" in line:
            if self.active:
                self.truncated_lines.append(self.block_line)
                self.incomplete = True
            self.seen = True
            self.active = True
            self.block_line = lineno
            self.free_energy = None
            self.energy_without_entropy = None
            self.energy_sigma0 = None
            return
        if not self.active:
            return
        match = re.search(r"\bfree\s+energy\s+TOTEN\s*=\s*(\S+)", line)
        if match:
            self.free_energy = match.group(1)
        match = re.search(r"\benergy\s+without\s+entropy\s*=\s*(\S+)", line)
        if match:
            self.energy_without_entropy = match.group(1)
        match = re.search(r"energy\(sigma->0\)\s*=\s*(\S+)", line)
        if match:
            self.energy_sigma0 = match.group(1)
        if self.free_energy is not None and self.energy_without_entropy is not None and self.energy_sigma0 is not None:
            self.active = False

    def complete(self) -> bool:
        return self.seen and not self.active and not self.incomplete

    def finish(self, issues: list[str]) -> None:
        for line in self.truncated_lines:
            _issue(issues, line, "truncated FREE ENERGIE block")
        if self.active:
            _issue(issues, self.block_line, "truncated FREE ENERGIE block")


class _FrameBuilder:
    def __init__(self, start_line: int) -> None:
        self.start_line = start_line
        self.cell_seen = False
        self.cell_pending = 0
        self.cell_rows: list[tuple[str, str, str]] = []
        self.positions_active = False
        self.positions_seen = False
        self.position_rows: list[tuple[str, str, str]] = []
        self.force_rows: list[tuple[str, str, str]] = []
        self.stress: tuple[str, ...] | None = None
        self.energy = _EnergyTracker()
        self.temperature: str | None = None

    def _finish_positions(self) -> None:
        self.positions_active = False

    def feed(self, line: str, lineno: int) -> None:
        stripped = line.strip()
        if "direct lattice vectors" in line:
            self.cell_seen = True
            self.cell_pending = 3
            self.cell_rows = []
            return
        if self.cell_pending:
            if not stripped:
                return
            tokens = stripped.split()
            if len(tokens) >= 3 and all(_number(token) for token in tokens[:3]):
                self.cell_rows.append(tuple(tokens[:3]))  # type: ignore[arg-type]
                self.cell_pending -= 1
                return
            self.cell_pending = -1

        if "POSITION" in line and "TOTAL-FORCE" in line:
            self.positions_active = True
            self.positions_seen = True
            self.position_rows = []
            return
        if self.positions_active:
            if not stripped or stripped.startswith("-"):
                return
            if stripped.lower().startswith("total drift"):
                self._finish_positions()
            else:
                tokens = stripped.split()
                if len(tokens) >= 6 and all(_number(token) for token in tokens[:6]):
                    self.position_rows.append((tokens[0], tokens[1], tokens[2]))
                    self.force_rows.append((tokens[3], tokens[4], tokens[5]))
                    return
                self._finish_positions()

        self.energy.feed(line, lineno)
        temperature = re.search(r"\btemperature\s+(\S+)\s*K\b", line)
        if temperature:
            self.temperature = temperature.group(1)

    def set_stress(self, stress: tuple[str, ...]) -> None:
        self.stress = stress

    def frame(self, index: int, expected_ions: int | None) -> OutcarFrame | None:
        if self.cell_pending or self.positions_active or not self.position_rows or not self.energy.complete():
            return None
        if self.cell_seen and len(self.cell_rows) != 3:
            return None
        if len(self.force_rows) != len(self.position_rows):
            return None
        if expected_ions is not None and len(self.position_rows) != expected_ions:
            return None
        return OutcarFrame(
            index,
            tuple(self.cell_rows) if self.cell_seen else None,
            tuple(self.position_rows),
            tuple(self.force_rows),
            self.stress,
            self.energy.free_energy,
            self.energy.energy_without_entropy,
            self.energy.energy_sigma0,
            self.temperature,
        )

    def count_issue(self, index: int, expected_ions: int | None) -> str | None:
        if not self.positions_seen:
            return None
        actual = len(self.position_rows)
        if expected_ions is not None and actual != expected_ions:
            return f"frame {index}: POSITION/TOTAL-FORCE row count {actual}, expected {expected_ions}"
        if self.positions_active:
            return f"frame {index}: incomplete POSITION/TOTAL-FORCE block after {actual} rows"
        return None

    def boundary_complete(self) -> bool:
        return self.stress is not None or (self.positions_seen and not self.positions_active)

    def has_data(self) -> bool:
        return self.cell_seen or self.positions_seen or self.stress is not None


class _ElasticTracker:
    def __init__(self) -> None:
        self.heading: str | None = None
        self.rows: list[tuple[str, ...]] = []

    def _finish(self, blocks: list[ElasticModuliBlock], issues: list[str], line: int) -> None:
        if self.heading is None:
            return
        if len(self.rows) == 6:
            blocks.append(ElasticModuliBlock(self.heading, tuple(self.rows)))
        else:
            _issue(issues, line, "truncated elastic-moduli block")
        self.heading = None
        self.rows = []

    def feed(self, line: str, lineno: int, blocks: list[ElasticModuliBlock], issues: list[str]) -> None:
        upper = line.upper()
        heading = next((item for item in _ELASTIC_HEADINGS if item in upper), None)
        if heading is not None:
            self._finish(blocks, issues, lineno)
            self.heading = line.rstrip("\r\n")
            return
        if self.heading is None:
            return
        tokens = line.split()
        if len(tokens) >= 7 and not _number(tokens[0]) and all(_number(token) for token in tokens[1:7]):
            self.rows.append(tuple(tokens[1:7]))
        elif len(tokens) >= 6 and all(_number(token) for token in tokens[:6]):
            self.rows.append(tuple(tokens[:6]))
        if len(self.rows) == 6:
            self._finish(blocks, issues, lineno)

    def finish(self, blocks: list[ElasticModuliBlock], issues: list[str], line: int) -> None:
        self._finish(blocks, issues, line)


def _is_separator(stripped: str) -> bool:
    return len(stripped) >= 4 and set(stripped) == {"-"}


class _MagnetizationTracker:
    """Track the final ``magnetization (x)`` block and its per-ion total moments.

    Only the collinear ``magnetization (x)`` table is parsed. Each new header
    resets the candidate so that the last block in the file always wins, and a
    truncated or malformed final block (from a killed job) yields ``None``. A
    ``magnetization (y)`` or ``(z)`` header seen after the final ``(x)`` header
    sets :attr:`noncollinear`, signalling that the ``(x)`` totals are only one
    projection of a noncollinear moment.
    """

    _IDLE = 0
    _AWAIT_ION_HEADER = 1
    _AWAIT_OPEN_SEP = 2
    _ROWS = 3
    _AWAIT_TOT = 4
    _DONE = 5

    def __init__(self) -> None:
        self.seen = False
        self.header_line = 0
        self._phase = self._IDLE
        self._rows: list[float] = []
        self._failed = False
        self.totals: tuple[float, ...] | None = None
        self.noncollinear = False

    def _reset(self, lineno: int) -> None:
        self.seen = True
        self.header_line = lineno
        self._phase = self._AWAIT_ION_HEADER
        self._rows = []
        self._failed = False
        self.totals = None
        self.noncollinear = False

    def feed(self, line: str, lineno: int) -> None:
        lowered = line.lower()
        if "magnetization (x)" in lowered:
            self._reset(lineno)
            return
        if self.seen and ("magnetization (y)" in lowered or "magnetization (z)" in lowered):
            self.noncollinear = True
            return
        if self._phase in (self._IDLE, self._DONE):
            return
        stripped = line.strip()
        if self._phase == self._AWAIT_ION_HEADER:
            if stripped.startswith("# of ion"):
                self._phase = self._AWAIT_OPEN_SEP
            return
        if self._phase == self._AWAIT_OPEN_SEP:
            if _is_separator(stripped):
                self._phase = self._ROWS
            elif stripped:
                self._failed = True
            return
        if self._phase == self._ROWS:
            if not stripped:
                return
            if _is_separator(stripped):
                self._phase = self._AWAIT_TOT
                return
            self._feed_row(stripped)
            return
        if self._phase == self._AWAIT_TOT:
            if not stripped:
                return
            if stripped.lower().startswith("tot") and not self._failed and self._rows:
                self.totals = tuple(self._rows)
            self._phase = self._DONE

    def _feed_row(self, stripped: str) -> None:
        tokens = stripped.split()
        if len(tokens) < 3 or not tokens[0].isdigit() or int(tokens[0]) <= 0:
            self._failed = True
            return
        if not all(_number(token) for token in tokens[1:]):
            self._failed = True
            return
        self._rows.append(float(tokens[-1]))

    def finish(self, issues: list[str]) -> None:
        if self.seen and self.totals is None:
            _issue(issues, self.header_line, "malformed magnetization (x) block")


def _iter_parsed_frames(
    lines: Iterator[str],
    *,
    issues: list[str] | None = None,
    stresses: list[tuple[str, ...]] | None = None,
    energy: _EnergyTracker | None = None,
    completion_evidence: list[str] | None = None,
    elastic: _ElasticTracker | None = None,
    elastic_blocks: list[ElasticModuliBlock] | None = None,
    magnetization: _MagnetizationTracker | None = None,
    line_counter: list[int] | None = None,
    nions: int | None = None,
) -> Iterator[OutcarFrame]:
    builder: _FrameBuilder | None = None
    next_index = 0
    expected_ions = nions
    last_line = 0
    elastic_issues = issues if issues is not None else []

    def finish_builder(current: _FrameBuilder) -> OutcarFrame | None:
        nonlocal next_index, expected_ions
        frame = current.frame(next_index, expected_ions)
        if frame is not None:
            if expected_ions is None:
                expected_ions = len(frame.positions or ())
            next_index += 1
            return frame
        if issues is not None and current.has_data():
            message = current.count_issue(next_index, expected_ions) or "truncated ionic frame"
            _issue(issues, current.start_line, message)
        return None

    for lineno, line in enumerate(lines, 1):
        last_line = lineno
        if line_counter is not None:
            line_counter[0] = lineno
        if energy is not None:
            energy.feed(line, lineno)
        if magnetization is not None:
            magnetization.feed(line, lineno)
        if completion_evidence is not None and _COMPLETION in line:
            completion_evidence.append(line.rstrip("\r\n"))
        stress = _stress_row(line)
        if stress is not None:
            if stresses is not None:
                stresses.append(stress)
            if builder is None:
                builder = _FrameBuilder(lineno)
            elif builder.boundary_complete():
                frame = finish_builder(builder)
                if frame is not None:
                    yield frame
                builder = _FrameBuilder(lineno)
            builder.set_stress(stress)
        if elastic is not None and elastic_blocks is not None:
            elastic.feed(line, lineno, elastic_blocks, elastic_issues)

        starts_section = "direct lattice vectors" in line or ("POSITION" in line and "TOTAL-FORCE" in line)
        if builder is not None and starts_section and builder.positions_seen and not builder.positions_active:
            frame = finish_builder(builder)
            if frame is not None:
                yield frame
            builder = _FrameBuilder(lineno)

        if _FRAME_MARKER.search(line):
            if builder is not None:
                frame = finish_builder(builder)
                if frame is not None:
                    yield frame
            builder = _FrameBuilder(lineno)
            continue
        if builder is not None:
            builder.feed(line, lineno)

    if builder is not None:
        frame = finish_builder(builder)
        if frame is not None:
            yield frame
    if elastic is not None and elastic_blocks is not None:
        elastic.finish(elastic_blocks, elastic_issues, last_line + 1)


class OutcarFile:
    """Lazy OUTCAR metadata reader whose scans reopen the source each time.

    OUTCAR paths, including compressed paths, are accepted by deliberate
    forward-streaming divergence from :class:`~httk.io.vasp.wavecar.WavecarFile`;
    random frame access re-streams the file. Construction validates only that
    the path exists. Prologue and full scans are lazy: the first prologue access
    scans to the ionic marker and can traverse the whole file when no marker
    exists. The full pass streams the source once and caches summary fields,
    all stress rows, and all elastic-moduli blocks. No source handle is
    retained, so :meth:`close` only marks this lazy object closed. The public
    :attr:`path` property returns the source filename.

    :param filename: Filesystem path to an OUTCAR, optionally compressed.
    """

    def __init__(self, filename: str | os.PathLike[str]) -> None:
        if not isinstance(filename, (str, os.PathLike)):
            raise TypeError("OUTCAR source must be a filesystem filename, not a live stream.")
        path = os.fsdecode(os.fspath(filename))
        if not os.path.isfile(path):
            raise FileNotFoundError(f"OUTCAR file does not exist: {path}")
        self._filename = path
        self._closed = False
        self._prologue: _Prologue | None = None
        self._full: _FullPass | None = None

    @property
    def closed(self) -> bool:
        """Whether this lazy reader has been closed."""
        return self._closed

    def close(self) -> None:
        """Mark this lazy object closed; it owns no persistent stream."""
        self._closed = True

    @property
    def path(self) -> str:
        """Return the source filename used to construct this lazy reader."""
        self._require_open()
        return self._filename

    def __enter__(self) -> Self:
        self._require_open()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def _require_open(self) -> None:
        if self.closed:
            raise ValueError("Cannot access a closed OUTCAR file.")

    def _ensure_prologue(self) -> _Prologue:
        if self._prologue is not None:
            return self._prologue

        version = ""
        nions: int | None = None
        ions_per_type: tuple[int, ...] | None = None
        parameters: dict[str, str] = {}
        parameter_lines: dict[str, int] = {}
        titles: list[str] = []
        issues: list[str] = []
        with source_lines(self._filename) as (lines, _raw):
            for lineno, line in enumerate(lines, 1):
                if _IONIC_MARKER.search(line):
                    break
                if not version:
                    version = next((token for token in line.split() if _VERSION_TOKEN.match(token)), "")
                nions_match = re.search(r"\bNIONS\s*=\s*(\d+)", line)
                if nions is None and nions_match is not None:
                    nions = int(nions_match.group(1))
                ions_match = re.search(r"\bions\s+per\s+type\s*=\s*(.*)$", line, re.IGNORECASE)
                if ions_per_type is None and ions_match is not None:
                    tokens = ions_match.group(1).split()
                    try:
                        ions_per_type = tuple(int(token) for token in tokens)
                    except ValueError:
                        _issue(issues, lineno, "unparseable ions per type")
                titel = _TITEL.match(line)
                if titel:
                    title = titel.group(1).strip()
                    if title not in titles:  # v1 deduplicated exact title strings, preserving first-seen order.
                        titles.append(title)
                for key in _PARAMETERS:
                    assignment = re.search(rf"\b{key}\s*=", line)
                    if assignment is None or key in parameters:
                        continue
                    value = re.search(rf"\b{key}\s*=\s*(\S+?)(?=\s|;|$)", line)
                    if value is None:
                        _issue(issues, lineno, f"unparseable {key} parameter")
                        continue
                    token = value.group(1)
                    parameters[key] = token
                    parameter_lines[key] = lineno
                    if key in _NUMERIC_PARAMETERS and re.fullmatch(_NUMBER, token) is None:
                        _issue(issues, lineno, f"unparseable {key} parameter value {token!r}")

        gga = parameters.get("GGA")
        xc_tag = parameters.get("LEXCH") if gga == "--" else gga
        xc = vasp_xc_tags.get(xc_tag) if xc_tag is not None else None
        if xc_tag is not None and xc is None:
            _issue(issues, parameter_lines.get("LEXCH" if gga == "--" else "GGA", 0), f"unknown XC tag {xc_tag!r}")

        self._prologue = _Prologue(
            version,
            _version_numbers(version),
            nions,
            ions_per_type,
            MappingProxyType(parameters),
            xc,
            tuple(titles),
            tuple(issues),
            MappingProxyType(parameter_lines),
        )
        return self._prologue

    def _ensure_full(self) -> _FullPass:
        if self._full is not None:
            return self._full
        prologue = self._ensure_prologue()
        issues = list(prologue.issues)
        energy = _EnergyTracker()
        completion_evidence: list[str] = []
        stresses: list[tuple[str, ...]] = []
        elastic_blocks: list[ElasticModuliBlock] = []
        elastic = _ElasticTracker()
        magnetization = _MagnetizationTracker()
        line_counter = [0]
        nframes = 0
        last_frame: OutcarFrame | None = None

        with source_lines(self._filename) as (lines, _raw):
            for frame in _iter_parsed_frames(
                iter(lines),
                issues=issues,
                stresses=stresses,
                energy=energy,
                completion_evidence=completion_evidence,
                elastic=elastic,
                elastic_blocks=elastic_blocks,
                magnetization=magnetization,
                line_counter=line_counter,
                nions=prologue.nions,
            ):
                nframes += 1
                last_frame = frame

        energy.finish(issues)
        magnetization.finish(issues)
        if not completion_evidence:
            _issue(issues, line_counter[0] + 1, "missing completion footer")
        self._full = _FullPass(
            FinalEnergies(
                energy.free_energy,
                energy.energy_without_entropy,
                energy.energy_sigma0,
                energy.seen,
            ),
            bool(completion_evidence),
            tuple(completion_evidence),
            tuple(issues),
            nframes,
            last_frame,
            tuple(stresses),
            tuple(elastic_blocks),
            magnetization.totals,
            magnetization.noncollinear,
        )
        return self._full

    @property
    def version_string(self) -> str:
        """Return the VASP version string found during the prologue scan."""
        self._require_open()
        return self._ensure_prologue().version_string

    @property
    def version_numbers(self) -> tuple[int, ...]:
        """Return the numeric components of the VASP version string."""
        self._require_open()
        return self._ensure_prologue().version_numbers

    @property
    def parameters(self) -> Mapping[str, str]:
        """Return the first recognized VASP parameter lexeme for each parameter."""
        self._require_open()
        return self._ensure_prologue().parameters

    @property
    def ions_per_type(self) -> tuple[int, ...] | None:
        """Return the number of ions for each potential type, when reported."""
        self._require_open()
        return self._ensure_prologue().ions_per_type

    @property
    def xc(self) -> str | None:
        """Return the recognized exchange-correlation description, when available."""
        self._require_open()
        return self._ensure_prologue().xc

    @property
    def potcar_titles(self) -> tuple[str, ...]:
        """Return distinct POTCAR titles in first-seen order."""
        self._require_open()
        return self._ensure_prologue().potcar_titles

    def frames(self) -> Iterator[OutcarFrame]:
        """Stream complete ionic frames without retaining the sequence.

        :return: An iterator yielding one :class:`OutcarFrame` at a time.
        """
        self._require_open()
        return self._frames()

    def _frames(self) -> Iterator[OutcarFrame]:
        self._require_open()
        prologue = self._ensure_prologue()
        with source_lines(self._filename) as (lines, _raw):
            yield from _iter_parsed_frames(iter(lines), nions=prologue.nions)

    def frame(self, index: int) -> OutcarFrame | None:
        """Return one frame by rescanning from the start of the file.

        :param index: Zero-based frame index.
        :return: The requested frame, or ``None`` when it is beyond the file.
        :raises ValueError: If ``index`` is negative or not an integer.
        """
        if not isinstance(index, int) or index < 0:
            raise ValueError("OUTCAR frame index must be a non-negative integer.")
        # ponytail: frame(i) is O(file); add a byte-offset index only if random access matters.
        return next(islice(self.frames(), index, index + 1), None)

    @property
    def final_energies(self) -> FinalEnergies:
        """Return energies from the last complete or partial energy block."""
        self._require_open()
        return self._ensure_full().final_energies

    @property
    def nframes(self) -> int:
        """Return the number of complete ionic frames after the full pass."""
        self._require_open()
        return self._ensure_full().nframes

    @property
    def last_frame(self) -> OutcarFrame | None:
        """Return the last complete ionic frame, when one exists."""
        self._require_open()
        return self._ensure_full().last_frame

    def stresses(self) -> tuple[tuple[str, ...], ...]:
        """Return all six-token ``in kB`` rows in file order.

        :return: Stress rows retaining their source numeric lexemes.
        """
        self._require_open()
        return self._ensure_full().stresses

    @property
    def elastic_moduli(self) -> tuple[ElasticModuliBlock, ...]:
        """Return parsed elastic-moduli tables in source order."""
        self._require_open()
        return self._ensure_full().elastic_moduli

    @property
    def magnetization(self) -> tuple[float, ...] | None:
        """Return per-ion total magnetic moments from the final ``magnetization (x)`` block.

        The values are the last (``tot``) column of each ion row in the last
        ``magnetization (x)`` block in the file, in Bohr magnetons. The
        per-orbital columns and the ``magnetization (y)`` / ``(z)`` blocks
        themselves are out of scope, so for a noncollinear run these values are
        only the *x* projection of each moment. A caller treating them as a
        collinear axis projection must check :attr:`noncollinear_magnetization`
        first. A malformed or truncated final block never raises: it yields
        ``None`` and records an entry in :attr:`issues`, so an OUTCAR from a
        killed job is ordinary input.

        :return: The per-ion total moments, or ``None`` for a non-spin-polarized
            run or an unusable final block.
        """
        self._require_open()
        return self._ensure_full().magnetization

    @property
    def noncollinear_magnetization(self) -> bool:
        """Whether a ``magnetization (y)`` or ``(z)`` block follows the final ``(x)`` block.

        When ``True`` the :attr:`magnetization` values are only the *x*
        projection of a noncollinear moment and must not be treated as a
        collinear axis magnitude.

        :return: ``True`` for a noncollinear final block, ``False`` otherwise,
            including when the file has no magnetization block.
        """
        self._require_open()
        return self._ensure_full().noncollinear_magnetization

    @property
    def completed(self) -> bool:
        """Whether the source contains completion-footer evidence."""
        self._require_open()
        return self._ensure_full().completed

    @property
    def completion_evidence(self) -> tuple[str, ...]:
        """Return completion-footer lines found during the full pass."""
        self._require_open()
        return self._ensure_full().completion_evidence

    @property
    def issues(self) -> tuple[str, ...]:
        """Return parsing issues collected during the available scans."""
        self._require_open()
        return self._ensure_full().issues


def read_outcar(source: Any) -> dict[str, Any]:
    """Return a lazy OUTCAR payload; only filesystem filenames are accepted.

    :param source: Filesystem path to an OUTCAR, optionally compressed.
    :return: A neutral payload containing the lazy OUTCAR reader.
    :raises TypeError: If ``source`` is not a filesystem path.
    :raises FileNotFoundError: If the path does not exist.
    """
    return {"format": "vasp-outcar", "outcar": OutcarFile(source)}
