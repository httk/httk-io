"""Header-only reader for VASP POTCAR and POTCAR.summary files.

Only potential header metadata is retained. The reader never returns or
otherwise exposes the full POTCAR text, whose body can contain licensed data.
``POTCAR.summary`` files made by concatenating those headers are supported.
"""

import re
from typing import Any

from ._text import source_lines


def _assignment(line: str, key: str) -> str | None:
    match = re.search(rf"\b{key}\s*=\s*([^;\s]+)", line)
    return match.group(1) if match else None


def _symbol(title: str) -> str:
    parts = title.split()
    return parts[1].split("_", 1)[0] if len(parts) > 1 else ""


def read_potcar_summary(source: Any) -> dict[str, Any]:
    """Extract one metadata mapping per ``TITEL`` header without retaining POTCAR text.

    Only header metadata is returned; the potential body is never retained or
    exposed because it can contain licensed data. Concatenated
    ``POTCAR.summary`` headers are accepted.

    :param source: POTCAR or POTCAR.summary filename, text stream, or iterable of source lines.
    :return: A neutral payload containing one metadata mapping per potential.
    """
    potentials: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None

    with source_lines(source) as (lines, _raw):
        for line in lines:
            titel_match = re.search(r"\bTITEL\s*=\s*(.*?)\s*$", line)
            if titel_match:
                if current is not None:
                    potentials.append(current)
                titel = titel_match.group(1).strip()
                current = {
                    "titel": titel,
                    "symbol": _symbol(titel),
                    "zval": "",
                    "pomass": "",
                    "enmax": "",
                    "lexch": None,
                }
                continue
            if current is None:
                continue
            current["zval"] = _assignment(line, "ZVAL") or current["zval"]
            current["pomass"] = _assignment(line, "POMASS") or current["pomass"]
            current["enmax"] = _assignment(line, "ENMAX") or current["enmax"]
            current["lexch"] = _assignment(line, "LEXCH") or current["lexch"]

    if current is not None:
        potentials.append(current)
    return {"format": "vasp-potcar", "potentials": potentials}
