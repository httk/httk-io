from httk.core.register import register_reader

register_reader(
    name="cif",
    reader="httk.io.cif:read_cif_asus",
    extensions=(".cif",),
)

from httk.core.register import register_writer

register_writer(
    name="cif",
    writer="httk.io.cif.cif_writer:_write_cif_payload",
    format="cif",
    extensions=(".cif",),
)

register_writer(
    name="poscar",
    writer="httk.io.vasp.poscar_writer:_write_poscar_payload",
    format="vasp-poscar",
    extensions=(".poscar", ".vasp"),
    filenames=("POSCAR", "CONTCAR"),
)

register_reader(
    name="mcif",
    reader="httk.io.cif:read_mcif_asus",
    extensions=(".mcif",),
)

register_reader(
    name="poscar",
    reader="httk.io.vasp:read_poscar",
    extensions=(".poscar", ".vasp"),
    filenames=("POSCAR", "CONTCAR"),
)
