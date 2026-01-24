from httk.core import register_loader

register_loader(
    name="cif",
    loader="httk.io.cif:read_cif",
    extensions=(".cif",),
)
