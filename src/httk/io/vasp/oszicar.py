"""String-preserving reader for VASP OSZICAR files.

Electronic iterations are attached to the following ionic summary, as VASP
writes them. If the file ends before that summary, a final entry with
``n``, ``F``, ``E0``, ``dE``, and ``mag`` all set to ``None`` preserves the
trailing electronic block.
"""

import re
from typing import Any

from ._text import source_lines

_COLON_LINE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_-]*)\s*:\s*(.*)$")
_IONIC_START = re.compile(r"^\s*(\d+)\b")


def _issue(issues: list[str], lineno: int, message: str) -> None:
    issues.append(f"line {lineno}: {message}")


def _value(line: str, key: str) -> str | None:
    match = re.search(rf"\b{key}\s*=\s*(\S+)", line)
    return match.group(1) if match else None


def _electronic(line: str) -> dict[str, Any] | None:
    match = _COLON_LINE.match(line)
    if match is None:
        return None
    scheme, rest = match.groups()
    fields = rest.split()
    if len(fields) not in (6, 7) or not fields[0].isdigit():
        return None
    return {
        "scheme": scheme,
        "n": int(fields[0]),
        "E": fields[1],
        "dE": fields[2],
        "d_eps": fields[3],
        "ncg": fields[4],
        "rms": fields[5],
        "rms_c": fields[6] if len(fields) == 7 else None,
    }


def read_oszicar(source: Any) -> dict[str, Any]:
    """Read OSZICAR text without converting numeric lexemes.

    Electronic iterations are attached to the following ionic summary. A
    trailing electronic block is retained as an incomplete final entry when no
    summary follows it.

    :param source: OSZICAR filename, text stream, or iterable of source lines.
    :return: A neutral payload containing ionic steps and parsing issues.
    """
    ionic_steps: list[dict[str, Any]] = []
    pending: list[dict[str, Any]] = []
    pending_start: int | None = None
    pending_last: int | None = None
    issues: list[str] = []

    with source_lines(source) as (lines, _raw):
        for lineno, line in enumerate(lines, 1):
            stripped = line.strip()
            if not stripped or (stripped.startswith("N ") and "d eps" in stripped):
                continue

            ionic_match = _IONIC_START.match(line)
            if ionic_match:
                n = int(ionic_match.group(1))
                boundary_pending = pending
                boundary_start = pending_start
                boundary_last = pending_last
                pending = []
                pending_start = None
                pending_last = None
                values = {key: _value(line, key) for key in ("F", "E0")}
                d_e = re.search(r"\bd\s*E\s*=\s*(\S+)", line)
                values["dE"] = _value(line, "dE") or (d_e.group(1) if d_e else None)
                # MD summaries commonly omit dE while retaining T/E/F/E0.
                if any(values[key] is None for key in ("F", "E0")):
                    if boundary_pending and boundary_start is not None and boundary_last is not None:
                        _issue(
                            issues,
                            lineno,
                            f"malformed ionic summary; dropped electronic iterations from lines "
                            f"{boundary_start}-{boundary_last}",
                        )
                    else:
                        _issue(issues, lineno, "malformed ionic summary")
                    continue
                ionic_steps.append(
                    {
                        "n": n,
                        "F": values["F"],
                        "E0": values["E0"],
                        "dE": values["dE"],
                        "mag": _value(line, "mag"),
                        "electronic": boundary_pending,
                    }
                )
                pending = []
                continue

            if _COLON_LINE.match(line):
                electronic = _electronic(line)
                if electronic is None:
                    _issue(issues, lineno, "malformed electronic iteration")
                else:
                    pending.append(electronic)
                    if pending_start is None:
                        pending_start = lineno
                    pending_last = lineno
                continue

            _issue(issues, lineno, "unrecognized line")

    if pending:
        ionic_steps.append({"n": None, "F": None, "E0": None, "dE": None, "mag": None, "electronic": pending})
    return {"format": "vasp-oszicar", "ionic_steps": ionic_steps, "issues": issues}
