"""Lazy directory-level access to a VASP calculation's registered files.

``VASPOutputs`` is deliberately not a registered reader: a directory is a
container of inputs, not one VASP file.  Payload-returning properties retain no
open handles; the lazy OUTCAR and XDATCAR objects are owned and closed here.
"""

import os
from pathlib import Path
from types import TracebackType
from typing import Any, Self

from httk.core.datastream import compression

from .oszicar import read_oszicar
from .outcar import OutcarFile
from .poscar_reader import read_poscar
from .potcar import read_potcar_summary
from .xdatcar import XdatcarFile

_MISSING = object()


def _compression_suffixes() -> tuple[str, ...]:
    # The core API exposes codec names, while this directory probe needs the
    # registered codec suffixes.  Reading the registry keeps custom codecs in scope.
    return tuple(ext for codec in compression._registry.values() for ext in codec.extensions)


class VASPOutputs:
    """Lazily resolve standard VASP files in a calculation directory.

    File lookup probes each registered compression suffix after the plain
    filename. Payload readers are created only when their properties are first
    accessed. The lazy OUTCAR and XDATCAR children are owned here and closed by
    :meth:`close`; payload mappings do not retain open handles.

    :param directory: Filesystem directory containing VASP output files.
    """

    def __init__(self, directory: str | os.PathLike[str]) -> None:
        if not isinstance(directory, (str, os.PathLike)):
            raise TypeError("VASPOutputs requires a filesystem directory")
        self._directory = Path(directory)
        if not self._directory.is_dir():
            raise NotADirectoryError(self._directory)
        self._closed = False
        self._cache: dict[str, Any] = {}

    @property
    def closed(self) -> bool:
        """Whether this directory view has been closed."""
        return self._closed

    def close(self) -> None:
        """Close owned lazy file objects; repeated calls are harmless."""
        if self._closed:
            return
        for name in ("outcar", "xdatcar"):
            value = self._cache.get(name)
            if value is not None and value is not _MISSING:
                value.close()
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
            raise ValueError("Cannot access closed VASP outputs.")

    def _find(self, name: str) -> Path | None:
        for suffix in ("", *_compression_suffixes()):
            path = self._directory / f"{name}{suffix}"
            if path.is_file():
                return path
        return None

    def _get(self, name: str, loader: Any, filename: str) -> Any:
        self._check_open()
        if name not in self._cache:
            path = self._find(filename)
            self._cache[name] = _MISSING if path is None else (loader(path) if loader else None)
        value = self._cache[name]
        return None if value is _MISSING else value

    @property
    def poscar(self) -> dict[str, Any] | None:
        """Return the lazily loaded POSCAR payload, or ``None`` when absent."""
        return self._get("poscar", read_poscar, "POSCAR")

    @property
    def contcar(self) -> dict[str, Any] | None:
        """Return the lazily loaded CONTCAR payload, or ``None`` when absent."""
        return self._get("contcar", read_poscar, "CONTCAR")

    @property
    def outcar(self) -> OutcarFile | None:
        """Return the owned lazy OUTCAR reader, or ``None`` when absent."""
        return self._get("outcar", OutcarFile, "OUTCAR")

    @property
    def xdatcar(self) -> XdatcarFile | None:
        """Return the owned lazy XDATCAR reader, or ``None`` when absent."""
        return self._get("xdatcar", XdatcarFile, "XDATCAR")

    @property
    def oszicar(self) -> dict[str, Any] | None:
        """Return the lazily loaded OSZICAR payload, or ``None`` when absent."""
        return self._get("oszicar", read_oszicar, "OSZICAR")

    @property
    def potcar(self) -> dict[str, Any] | None:
        """Return the lazily loaded POTCAR summary, or ``None`` when absent."""
        self._check_open()
        if "potcar" not in self._cache:
            path = self._find("POTCAR") or self._find("POTCAR.summary")
            self._cache["potcar"] = _MISSING if path is None else read_potcar_summary(path)
        value = self._cache["potcar"]
        if value is _MISSING:
            return None
        if not isinstance(value, dict):
            raise TypeError("cached POTCAR payload is not a mapping")
        return value
