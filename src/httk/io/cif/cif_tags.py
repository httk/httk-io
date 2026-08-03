"""CIF data names after :func:`httk.io.cif.cif_reader.read_cif` has normalized them."""

from typing import TypedDict


class CifTags(TypedDict):
    """The normalized spellings of the special-purpose CIF data names this package consults."""

    structural_q: tuple[str, ...]
    structural_displacement_label: str
    structural_occupancy_label: str
    magnetic_cartesian_moment: tuple[str, ...]
    magnetic_fourier_coeff: str
    magnetic_fourier_label: str
    magnetic_ssg_name: str


def _tag(name: str) -> str:
    return name.lstrip("_").lower()


CIF_TAGS: CifTags = {
    "structural_q": tuple(_tag(f"_cell_wave_vector_{axis}") for axis in "xyz"),
    "structural_displacement_label": _tag("_atom_site_displace_Fourier.atom_site_label"),
    "structural_occupancy_label": _tag("_atom_site_occupancy_Fourier.atom_site_label"),
    "magnetic_cartesian_moment": tuple(_tag(f"_atom_site_moment.Cartn_{axis}") for axis in "xyz"),
    "magnetic_fourier_coeff": _tag("_atom_site_Fourier_wave_vector.q{}_coeff"),
    "magnetic_fourier_label": _tag("_atom_site_moment_Fourier.atom_site_label"),
    "magnetic_ssg_name": _tag("_space_group.magn_ssg_name"),
}
