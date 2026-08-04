"""Magnetic-CIF fixtures use CIF2 lists only; tables and triple-quoted strings remain unsupported."""

from fractions import Fraction
from pathlib import Path

import pytest

from httk.io.cif.cif_parser import parse_structural_modulation
from httk.io.cif.cif_reader import read_cif
from httk.io.cif.mcif_parser import (
    _parse_modulation,
    _parse_xyzt_op,
    cifblock_to_mag_asu,
    extract_fourier,
    mag_asus_from_mcif_file,
    single_mag_asu_from_mcif_file,
)

FIXTURES = Path(__file__).with_name("fixtures")


def test_centered_magnetic_group_composes_translation_and_time_reversal() -> None:
    data = single_mag_asu_from_mcif_file(FIXTURES / "magnetic_centered.mcif")

    assert data["format"] == "mcif"
    assert data["moment_basis"] == "crystalaxis"
    assert data["magmoms_exact"] == (("1", "0", "0"),)
    assert data["centerings_xyz"] == ("x,y,z,+1", "x+1/2,y+1/2,z+1/2,-1")
    assert "spin_basis" not in data
    assert "magmoms" not in data
    assert data["cell_parameters_exact"] == ("1",) * 3 + ("90",) * 3
    assert data["positions_exact"] == [("0", "0", "0")]
    assert data["symops_xyz"] == ("x,y,z,+1",)
    assert "basis" not in data
    assert "positions" not in data
    assert "symops" not in data


def test_algebraic_ssg_operations_use_structural_q_and_mark_magnetic_modulation() -> None:
    data = single_mag_asu_from_mcif_file(FIXTURES / "magnetic_ssg.mcif")

    assert data["format"] == "mcif"
    assert data["moment_basis"] == "crystalaxis"
    assert data["magmoms_exact"] == (("0", "0", "1"),)
    assert data["centerings_xyz"] == ("x1,x2,x3,x4,+1",)
    assert data["symops_xyz"] == ("x1,x2,x3,x4,+1",)
    assert "symops" not in data
    assert data["incomm"]["magnetic_q"] == [[0.0, 0.0, pytest.approx(1 / 3)]]
    assert data["incomm"]["has_magnetic_modulation"] is True
    assert data["incomm"]["magnetic_modulated_atoms"] == ["Fe1"]
    assert data["space_group_name_hm"] is None


def test_cartesian_moments_use_cartn_columns_without_crystal_axis_conversion() -> None:
    data = single_mag_asu_from_mcif_file(FIXTURES / "magnetic_cartesian.mcif")

    # A 2 x 3 x 4 orthogonal cell makes the intended Cartesian result hand-checkable.
    assert data["format"] == "mcif"
    assert data["moment_basis"] == "cartesian"
    assert data["magmoms_exact"] == (("1", "2", "3"),)
    assert data["centerings_xyz"] == ("x,y,z,+1",)


def test_cartesian_moments_are_not_converted_through_a_hexagonal_crystal_basis() -> None:
    data = single_mag_asu_from_mcif_file(FIXTURES / "magnetic_cartesian_hexagonal.mcif")

    # With gamma=120°, an incorrect crystal-axis conversion of (1, 2, 3) would give
    # (0, sqrt(3), 3). Cartn values are already Cartesian and must remain unchanged.
    assert data["format"] == "mcif"
    assert data["moment_basis"] == "cartesian"
    assert data["magmoms_exact"] == (("1", "2", "3"),)
    assert data["centerings_xyz"] == ("x,y,z,+1",)
    assert data["positions_exact"] == [("0", "0", "0")]
    assert data["symops_xyz"] == ("x,y,z,+1",)


def test_private_xyzt_parser_rejects_invalid_time_reversal() -> None:
    with pytest.raises(ValueError, match="Invalid time-reversal"):
        _parse_xyzt_op("1x,1y,1z,+3", use_fractions=True, time_reversal_convention="spglib")


def test_k_vector_and_fourier_coefficients_are_exact() -> None:
    _name, block = read_cif(FIXTURES / "magnetic_kvector.mcif", allow_cif2=True)[0][0]

    _structural_q, magnetic_q, *_rest = _parse_modulation(block)
    assert magnetic_q == [[Fraction(1, 8), Fraction(0), Fraction(1, 3)]]
    assert extract_fourier(block) == ([(Fraction(1, 8), Fraction(0), Fraction(1, 3))], [(1,)])


def test_structural_modulation_tags_are_detected_with_normalized_names() -> None:
    _name, block = read_cif(FIXTURES / "structural_modulation.cif")[0][0]

    q_vectors, dimension, has_modulation, labels = parse_structural_modulation(block)
    assert q_vectors == [[Fraction(0), Fraction(0), Fraction(1, 3)]]
    assert dimension == 1
    assert has_modulation is True
    assert labels == ["Fe1", "Fe2"]


def test_short_loop_names_the_loop_and_column_counts(tmp_path: Path) -> None:
    malformed = tmp_path / "short.cif"
    malformed.write_text("data_bad\nloop_\n_a\n_b\none\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"loop.*2.*1"):
        read_cif(malformed)


def test_moment_column_mismatch_can_be_rejected_or_treated_as_nonmagnetic(tmp_path: Path) -> None:
    _name, malformed = read_cif(FIXTURES / "magnetic_cartesian.mcif")[0][0]
    malformed["atom_site_moment.cartn_z"] = []

    assert cifblock_to_mag_asu(malformed)["magmoms_exact"] is None
    with pytest.raises(ValueError, match="moment columns"):
        cifblock_to_mag_asu(malformed, error_on_nonmag=True)


def test_moments_fill_unlisted_atom_sites_with_exact_zero_tokens() -> None:
    _name, block = read_cif(FIXTURES / "magnetic_centered.mcif")[0][0]
    block["atom_site_label"].append("Fe2")
    block["atom_site_type_symbol"].append("Fe")
    block["atom_site_fract_x"].append("1/2")
    block["atom_site_fract_y"].append("1/2")
    block["atom_site_fract_z"].append("1/2")

    data = cifblock_to_mag_asu(block)

    assert data["labels"] == ["Fe1", "Fe2"]
    assert data["magmoms_exact"] == (("1", "0", "0"), ("0", "0", "0"))


def test_blocks_without_a_moment_loop_have_no_moment_payload() -> None:
    _name, block = read_cif(FIXTURES / "magnetic_centered.mcif")[0][0]
    for tag in (
        "atom_site_moment.label",
        "atom_site_moment.crystalaxis_x",
        "atom_site_moment.crystalaxis_y",
        "atom_site_moment.crystalaxis_z",
    ):
        block.pop(tag)

    data = cifblock_to_mag_asu(block)

    assert data["moment_basis"] is None
    assert data["magmoms_exact"] is None


def test_cif2_nested_and_empty_lists_reach_the_mcif_loader() -> None:
    data = mag_asus_from_mcif_file(FIXTURES / "cif2_lists.mcif")
    _name, block = read_cif(FIXTURES / "cif2_lists.mcif", allow_cif2=True)[0][0]

    assert len(data) == 1
    assert block["audit_creation_method"] == ["alpha", ["beta"], []]
