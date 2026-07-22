"""Smoke tests: the package and its handler registration package import cleanly."""


def test_import_httk_io():
    import httk.io

    assert hasattr(httk.io, "read_cif")


def test_import_cif_subpackage():
    import httk.io.cif as cif

    for name in cif.__all__:
        assert hasattr(cif, name), name


def test_import_handlers_package():
    # Importing the handler package must register the cif loader as a side effect.
    import httk.handlers.io  # noqa: F401
