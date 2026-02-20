from httk.core.load import register_loader

register_loader(
    name="cif",
    loader="httk.io.cif:read_cif",
    extensions=(".cif",),
)
