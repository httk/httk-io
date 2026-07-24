from httk.core.register import register_loader

register_loader(
    name="cif",
    loader="httk.io.cif:read_cif",
    extensions=(".cif",),
)

register_loader(
    name="poscar",
    loader="httk.io.vasp:read_poscar",
    extensions=(".poscar", ".vasp"),
    filenames=("POSCAR", "CONTCAR"),
)
