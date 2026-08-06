"""Small source adapters shared by VASP text readers."""

import os
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from httk.core import BytestreamFileView, TextstreamFileView


@contextmanager
def source_lines(
    source: Any, *, preserve_path: bool = False, capture_stream: bool = False
) -> Iterator[tuple[Iterable[str], str | None]]:
    """Yield source lines and, when available, the original decoded text."""
    opened: BytestreamFileView | TextstreamFileView | None = None
    if isinstance(source, (str, os.PathLike)):
        if preserve_path:
            opened = BytestreamFileView(Path(source))
            text = opened.read().decode("utf-8")
            lines: Iterable[str] = text.splitlines(keepends=True)
            raw: str | None = text
        else:
            opened = TextstreamFileView(Path(source))
            lines = opened
            raw = None
    elif capture_stream and hasattr(source, "read"):
        text = source.read()
        if isinstance(text, bytes):
            text = text.decode("utf-8")
        lines = text.splitlines(keepends=True)
        raw = text
    elif hasattr(source, "readline"):
        lines = source
        raw = None
    elif hasattr(source, "read"):
        text = source.read()
        if isinstance(text, bytes):
            text = text.decode("utf-8")
        lines = text.splitlines(keepends=True)
        raw = None
    else:
        lines = source
        raw = None
    try:
        yield lines, raw
    finally:
        if opened is not None:
            opened.close()
