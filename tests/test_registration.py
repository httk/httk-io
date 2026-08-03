"""The cif loader is discoverable through httk-core's registry."""

import httk.core


def test_cif_extension_registered():
    # Discovery runs on `import httk.core` and imports httk.registry.io,
    # which registers the ".cif" loader.
    assert ".cif" in httk.core.register.known_extensions()


def test_cif_loader_points_at_the_asu_reader():
    """``load`` yields interpreted asymmetric units, not the raw token tree.

    The low-level ``read_cif`` tokenizer is still exported for callers who want the tags
    verbatim; it is just not what the registry dispatches to, so that a ``.cif`` behaves
    like a ``POSCAR`` and can be handed straight to the structure builders.
    """
    spec = httk.core.register.loaders.require(".cif")
    assert spec.name == "cif"
    assert spec.handler == "httk.io.cif:read_cif_asus"


def test_cif_loader_resolves_to_callable():
    from httk.core._plugins import resolve_callable

    spec = httk.core.register.loaders.require(".cif")
    fn = resolve_callable(spec.handler)
    from httk.io.cif.cif_parser import read_cif_asus

    assert fn is read_cif_asus


def test_mcif_extension_registered_with_neutral_loader():
    spec = httk.core.register.loaders.require(".mcif")
    assert spec.name == "mcif"
    assert spec.handler == "httk.io.cif:read_mcif_asus"
