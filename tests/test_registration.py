"""The cif loader is discoverable through httk-core's registry."""

import httk.core


def test_cif_extension_registered():
    # Discovery runs on `import httk.core` and imports httk.handlers.io,
    # which registers the ".cif" loader.
    assert ".cif" in httk.core.register.known_extensions()


def test_cif_loader_points_at_read_cif():
    spec = httk.core.register.loaders.require(".cif")
    assert spec.name == "cif"
    assert spec.handler == "httk.io.cif:read_cif"


def test_cif_loader_resolves_to_callable():
    from httk.core._plugins import resolve_callable

    spec = httk.core.register.loaders.require(".cif")
    fn = resolve_callable(spec.handler)
    from httk.io.cif.cif_reader import read_cif

    assert fn is read_cif
