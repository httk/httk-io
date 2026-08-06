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

register_reader(
    name="oszicar",
    reader="httk.io.vasp:read_oszicar",
    extensions=(".oszicar",),
    filenames=("OSZICAR",),
)

register_reader(
    name="outcar",
    reader="httk.io.vasp:read_outcar",
    extensions=(".outcar",),
    filenames=("OUTCAR",),
)

register_reader(
    name="potcar",
    reader="httk.io.vasp:read_potcar_summary",
    extensions=(".potcar",),
    filenames=("POTCAR", "POTCAR.summary"),
)

register_reader(
    name="xdatcar",
    reader="httk.io.vasp:read_xdatcar",
    extensions=(".xdatcar",),
    filenames=("XDATCAR",),
)

register_reader(
    name="wavecar",
    reader="httk.io.vasp.wavecar:read_wavecar",
    extensions=(".wavecar",),
    filenames=("WAVECAR",),
)

register_reader(
    name="trajectory-jsonl",
    reader="httk.io.optimade_jsonl:read_trajectory_jsonl",
    extensions=(".jsonl",),
)

register_writer(
    name="wavecar",
    writer="httk.io.vasp.wavecar:_write_wavecar_payload",
    format="vasp-wavecar",
    extensions=(".wavecar",),
    filenames=("WAVECAR",),
)

register_writer(
    name="trajectory-jsonl",
    writer="httk.io.optimade_jsonl:_write_trajectory_jsonl_payload",
    format="httk-trajectory-jsonl",
    extensions=(".jsonl",),
)
