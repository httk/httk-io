#
#    The high-throughput toolkit (httk)
#    Copyright (C) 2012-2025 The httk AUTHORS
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import re
from collections.abc import Iterable
from fractions import Fraction
from typing import Any, cast

from httk.core import combined_precision

from ._xyz_expr import _parse_linear_expr, _parse_linear_expr_algebraic
from .cif_parser import (
    _cell_parameter_tokens,
    cif_exact_token,
    parse_asu_cell,
    parse_cif_float,
    parse_cif_fraction,
    parse_structural_modulation,
)
from .cif_reader import read_cif
from .cif_tags import CIF_TAGS


def extract_parent_q_basis(cifblock: dict[str, Any]) -> list[tuple[Fraction, Fraction, Fraction]] | None:
    """
    Return the parent propagation basis as a list of ``(kx, ky, kz)`` tuples,
    or ``None`` if it is not present.
    """
    k_vectors = cifblock.get('parent_propagation_vector.kxkykz')
    if not k_vectors:
        return None
    basis = []
    for row in k_vectors:
        # each row like ('0', '0', '1/3') or [0,0,0.333...]
        vector = tuple(parse_cif_fraction(v) for v in row)
        if len(vector) != 3 or any(value is None for value in vector):
            raise ValueError(f"Invalid parent propagation vector: {row!r}")
        basis.append(cast(tuple[Fraction, Fraction, Fraction], vector))
    return basis


def extract_fourier_coeffs(cifblock: dict[str, Any]) -> tuple[list[tuple[Any, ...]], int]:
    """
    Return ``(coeff_rows, m)`` where ``coeff_rows`` is a list of coefficient tuples
    ``(c1, c2, ..., cm)`` and ``m`` is the number of q-vectors detected (``>= 0``).

    All present ``q{i}_coeff`` columns are found, zipped row-wise with missing entries
    filled by zeros, and duplicate coefficient tuples are removed.
    """
    # discover which q*_coeff columns exist
    present_cols = []
    pattern = re.compile("^" + re.escape(CIF_TAGS['magnetic_fourier_coeff']).replace(r"\{\}", r"(\d+)") + "$")
    for key, col in cifblock.items():
        match = pattern.match(key)
        if match and col is not None:
            present_cols.append((int(match.group(1)), key, col))

    if not present_cols:
        return [], 0

    # All present columns must have same length (number of rows). We’ll be permissive: pad shorter ones with zeros.
    max_len = max(len(col) for (_, _, col) in present_cols)
    m = max(i for (i, _, _) in present_cols)

    # Build a dense matrix of size (max_len x m), filling missing cols with zeros
    # 1-based indexing externally; 0-based in list
    columns: list[list[Any] | None] = [None] * m
    for i, key, col in present_cols:
        # normalize to numeric (ints preferred) but accept rational/float
        def norm(x):
            # msCIF usually stores integer coeffs; permit '1', '-1', '0', '2/3', '0.5'
            s = str(x).strip()
            if "/" in s:
                return Fraction(s)
            try:
                v = int(s)
                return v
            except Exception:
                try:
                    return float(s)
                except Exception:
                    return s

        col_norm = [norm(v) for v in col]
        if len(col_norm) < max_len:
            col_norm = col_norm + [0] * (max_len - len(col_norm))
        columns[i - 1] = col_norm

    # Any missing columns among 1..m become zeros
    for idx in range(m):
        if columns[idx] is None:
            columns[idx] = [0] * max_len

    # transpose to rows
    rows = list(zip(*(cast(list[Any], column) for column in columns)))

    # deduplicate coefficient tuples
    coeff_rows = []
    seen = set()
    for r in rows:
        tup = tuple(r)
        # For dedup, coerce Fractions to a canonical string (since Fraction is hashable, this is optional)
        if tup not in seen:
            seen.add(tup)
            coeff_rows.append(tup)

    return coeff_rows, m


def extract_fourier(cifblock):
    """
    Return the ``(basis, coeffs)`` descriptor, or ``None`` if there is insufficient data.

    ``basis`` comes from ``parent_propagation_vector.kxkykz`` and ``coeffs`` are the
    unique coefficient tuples from ``atom_site_Fourier_wave_vector.q*_coeff``.
    """
    basis = extract_parent_q_basis(cifblock)
    coeff_rows, m = extract_fourier_coeffs(cifblock)

    if not basis or not coeff_rows:
        # Not enough to build a fourier descriptor
        return None

    if len(basis) < m:
        # If fewer basis vectors than coeff columns, pad missing q’s with (0,0,0)
        # (harmless for commensurability; or you can raise if you prefer strictness)
        basis = list(basis) + [(Fraction(0), Fraction(0), Fraction(0))] * (m - len(basis))
    elif len(basis) > m:
        # Truncate extra basis vectors (common if multiple parent k’s present but only q1..qm used)
        basis = list(basis[:m])

    return (basis, coeff_rows)


def _collect_k_from_fourier(fourier):
    """
    Collect the magnetic propagation k-vectors from a fourier descriptor.

    `fourier` is the ``(basis, coeff_rows)`` tuple returned by
    :func:`extract_fourier`; the basis vectors are the parent propagation
    k-vectors, returned here as a list of ``[kx, ky, kz]`` lists.
    """
    basis, _coeff_rows = fourier
    return [list(k) for k in basis]


def _parse_xyzt_op(op, use_fractions=False, time_reversal_convention="mcif"):
    """
    op: e.g. 'x-y,x,-z+1/2,-1'
    Returns (R, t, time) where R is 3x3 (list of rows), t is length-3 float list, time is +1/-1.
    """
    parts = [p.strip() for p in op.split(",")]
    if len(parts) != 4:
        raise ValueError(f"Unexpected op format: {op}")
    px, py, pz, ts = parts
    rx, tx = _parse_linear_expr(px, use_fractions=use_fractions)
    ry, ty = _parse_linear_expr(py, use_fractions=use_fractions)
    rz, tz = _parse_linear_expr(pz, use_fractions=use_fractions)
    try:
        time_value = int(ts)
    except ValueError as error:
        raise ValueError(f"Invalid time-reversal flag at end of operation: '{ts}' in {op}") from error
    if time_value not in (-1, 1):
        raise ValueError(f"Invalid time-reversal flag at end of operation: '{ts}' in {op}")
    if time_reversal_convention == "mcif":
        ts = time_value
    elif time_reversal_convention == "spglib":
        ts = int((1 - time_value) / 2)
    else:
        raise ValueError(f"Unrecognized time-reversal convention: {time_reversal_convention!r}")
    return (rx, ry, rz), (tx, ty, tz), ts


def _xyzt_symops_to_matrix(symops_xyz, use_fractions=False, time_reversal_convention="mcif"):
    return [_parse_xyzt_op(s, use_fractions, time_reversal_convention=time_reversal_convention) for s in symops_xyz]


def _compose_ops_with_centerings(ops, centerings):
    """
    ops:         list of (R, t, time_flag) from _space_group_symop_magn_operation
    centerings:  list of (c, time_c) where c is 3-vector fractional translation,
                 time_c is 0 or 1 (0 = no time reversal)

    Returns a new list of (R, t', time') where t' = t + c and time' = (time_flag + time_c)%2.
    """
    composed = []
    for R, t, time_flag in ops:
        for Rc, c, time_c in centerings:
            if Rc != ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                raise ValueError("Magnetic centering symmetry operation must not include rotation")
            t_new = (t[0] + c[0], t[1] + c[1], t[2] + c[2])
            time_new = (time_flag + time_c) % 2  # time_flag * time_c for -1/+1 convention
            composed.append((R, t_new, time_new))
    return composed


def _parse_moments(block, *, resolution=True) -> tuple[Any, ...] | None:
    """
    Extract magnetic moment tokens and their frame from an mcif block.

    Parameters
    ----------
    block: mcif block
    resolution : bool
        If True, return the combined precision of the moment tokens.

    Returns
    -------
    If resolution=False:
        moments, labels, moment_basis
    If resolution=True:
        moments, labels, moment_basis, magmom_precision

    ``magmom_precision`` is in μB for either frame: crystalaxis components are
    specified along the unit lattice axes.
    """

    def _get(name):
        v = block.get(name)
        return [] if v is None else list(v)

    def _len_ok(xs, ys, zs, n):
        return len(xs) == len(ys) == len(zs) == n and n > 0

    if 'atom_site_moment.label' not in block:
        return None

    labels = _get('atom_site_moment.label')
    n = len(labels)
    if n == 0:
        return None

    # Try crystal basis first
    xs = _get('atom_site_moment.crystalaxis_x')
    ys = _get('atom_site_moment.crystalaxis_y')
    zs = _get('atom_site_moment.crystalaxis_z')
    moment_basis = "crystalaxis"

    if not _len_ok(xs, ys, zs, n):
        xs, ys, zs = (_get(tag) for tag in CIF_TAGS['magnetic_cartesian_moment'])
        moment_basis = "cartesian"
        if not _len_ok(xs, ys, zs, n):
            return None if not resolution else (None, None, None, None)

    # Labels must be unique
    if n != len(set(labels)):
        raise ValueError("Non-equivalent sites share the same moment label in CIF data.")

    moments = []
    component_claims: list[object] = []

    for i in range(n):
        if resolution:
            _mx, mx_meta = parse_cif_float(xs[i], meta=True)
            _my, my_meta = parse_cif_float(ys[i], meta=True)
            _mz, mz_meta = parse_cif_float(zs[i], meta=True)

            component_claims.extend(
                [
                    mx_meta['precision'],
                    mx_meta['esd'],
                    my_meta['precision'],
                    my_meta['esd'],
                    mz_meta['precision'],
                    mz_meta['esd'],
                ]
            )

        moments.append((cif_exact_token(xs[i]), cif_exact_token(ys[i]), cif_exact_token(zs[i])))

    # Fast path: no grid resolution requested
    if not resolution:
        return moments, labels, moment_basis

    # Crystalaxis components are μB along the unit lattice axes, so this precision is in
    # μB for both moment bases.
    mag_res = combined_precision(component_claims)

    return moments, labels, moment_basis, mag_res


def _parse_alg_op(
    op: str, use_fractions: bool = False, time_reversal_convention: str = "mcif"
) -> tuple[Any, tuple[Any, Any, Any], int]:
    """
    Parse an msCIF `_space_group_symop_magn_ssg_operation.algebraic` string.

    Examples of `op`:
      'x1,x2,x3,x4,x5,x6,+1'
      'x1,x2,x3,x4+1/2,-1'
      'x1+1/2,x2+1/2,-x3,-x4,-1'
      'x1-x2,x1,x3+1/3,x4-1/6,x5,+1'

    Returns (R, t, time):
      R: 3x3 list of rows with entries in {-1,0,1}
      t: length-3 list of floats (or Fractions if use_fractions=True)
      time: +1 / -1 (or mapped for spglib)
    """
    parts = [p.strip() for p in op.split(",")]
    if len(parts) < 4:
        raise ValueError(f"Unexpected op format (need at least 3 coords + time): {op}")

    # last item is the time-reversal flag ('+1' / '-1')
    ts_str = parts[-1]
    coord_parts = parts[:-1]

    if len(coord_parts) < 3:
        raise ValueError(f"Need at least 3 coordinate expressions before time flag: {op}")

    # Only the first three coordinates define the 3D spatial mapping we return.
    px, py, pz = coord_parts[0], coord_parts[1], coord_parts[2]

    # Parse rows; only x1,x2,x3 are permitted to appear in the first three expressions.
    allowed = ('x1', 'x2', 'x3')
    rx, tx = _parse_linear_expr_algebraic(px, allowed_vars=allowed, use_fractions=use_fractions)
    ry, ty = _parse_linear_expr_algebraic(py, allowed_vars=allowed, use_fractions=use_fractions)
    rz, tz = _parse_linear_expr_algebraic(pz, allowed_vars=allowed, use_fractions=use_fractions)

    # Time reversal handling
    try:
        ts_val = int(ts_str)
        if ts_val not in (-1, 1):
            raise ValueError
    except Exception:
        raise ValueError(f"Invalid time-reversal flag at end of operation: '{ts_str}' in {op}")

    if time_reversal_convention == "mcif":
        time = ts_val
    elif time_reversal_convention == "spglib":
        # Map +1 -> 0 (no time-reversal), -1 -> 1 (time-reversal)
        time = int((1 - ts_val) // 2)
    else:
        raise ValueError("Unrecognized time reversal convention. Use 'mcif' or 'spglib'.")

    # Prepare translations as requested type
    if use_fractions:
        t = (tx, ty, tz)
    else:
        t = (float(tx), float(ty), float(tz))

    R = (rx, ry, rz)
    return R, t, time


def _alg_symops_to_matrix(symops_alg, use_fractions=False, time_reversal_convention="mcif"):
    return [_parse_alg_op(s, use_fractions, time_reversal_convention=time_reversal_convention) for s in symops_alg]


def _parse_mag_asu_cell(cifblock: dict[str, Any]) -> tuple[Any, ...]:
    asu = parse_asu_cell(cifblock)
    moments_result = _parse_moments(cifblock, resolution=True)
    if moments_result is None:
        return asu, None, None, None

    cif_moments, momlabels, moment_basis, magres = moments_result

    if cif_moments is None:
        return (
            asu,
            None,
            None,
            magres,
        )

    moments_map = {label: (mom[0], mom[1], mom[2]) for label, mom in zip(momlabels, cif_moments)}
    magmoms_exact = tuple(moments_map.get(i, ("0", "0", "0")) for i in asu.labels)

    return asu, magmoms_exact, moment_basis, magres


def _get_magnetic_fourier_info(cifblock: dict[str, Any]) -> tuple[bool, list[str]]:
    has = False
    atoms = set()

    labels = cifblock.get(CIF_TAGS['magnetic_fourier_label'])
    if labels:
        has = True
        atoms.update(labels)

    return has, sorted(atoms)


def _parse_modulation(cifblock: dict[str, Any]) -> tuple[Any, ...]:
    structural_q, mod_dim, has_struct_mod, struct_mod_atoms = parse_structural_modulation(cifblock)

    # Magnetic q
    magnetic_q = None

    # (A) magnetic superspace -> uses same q as structural superspace
    if cifblock.get(CIF_TAGS['magnetic_ssg_name']) is not None and structural_q:
        magnetic_q = structural_q

    # (B) commensurate magnetic propagation vector
    elif cifblock.get('parent_propagation_vector.kxkykz'):
        rows = cifblock['parent_propagation_vector.kxkykz']
        vectors: list[list[Fraction]] = []
        for row in rows:
            vector = [parse_cif_fraction(v) for v in row]
            if len(vector) != 3 or any(value is None for value in vector):
                raise ValueError(f"Invalid parent propagation vector: {row!r}")
            vectors.append(cast(list[Fraction], vector))
        magnetic_q = vectors

    # (C) Fourier-defined magnetic propagation vector
    else:
        fourier = extract_fourier(cifblock)
        if fourier:
            magnetic_q = _collect_k_from_fourier(fourier)

    has_mag_mod, mag_mod_atoms = _get_magnetic_fourier_info(cifblock)

    return structural_q, magnetic_q, mod_dim, has_struct_mod, has_mag_mod, struct_mod_atoms, mag_mod_atoms


def cifblock_to_mag_asu(cifblock: dict[str, Any], *, error_on_nonmag: bool = False) -> dict[str, Any]:
    (
        asu,
        magmoms_exact,
        moment_basis,
        magres,
    ) = _parse_mag_asu_cell(cifblock)
    structural_q, magnetic_q, mod_dim, has_struct_mod, has_mag_mod, struct_mod_atoms, mag_mod_atoms = _parse_modulation(
        cifblock
    )

    if error_on_nonmag and magmoms_exact is None:
        raise ValueError("Magnetic moment columns are missing or have mismatched lengths")

    # Exact numeric CIF tokens are Fractions, so their commensurability never relies on
    # a floating-point denominator guess.
    incomm = None
    if mod_dim > 0:
        incomm = {
            'structural_q': structural_q,
            'magnetic_q': magnetic_q,
            'mod_dim': mod_dim,
            'has_structural_modulation': has_struct_mod,
            'structural_modulated_atoms': struct_mod_atoms,
            'has_magnetic_modulation': has_mag_mod,
            'magnetic_modulated_atoms': mag_mod_atoms,
        }

    base_symops_xyz = cifblock.get('space_group_symop_magn_operation.xyz')
    if base_symops_xyz is None:
        base_symops_alg = cifblock.get('space_group_symop_magn_ssg_operation.algebraic')
        if base_symops_alg is None:
            raise ValueError("mcif block has no symmetry operations")
        raw_symops = cast(list[str], base_symops_alg)
        base_symops = _alg_symops_to_matrix(base_symops_alg, use_fractions=True, time_reversal_convention="spglib")
    else:
        raw_symops = cast(list[str], base_symops_xyz)
        base_symops = _xyzt_symops_to_matrix(base_symops_xyz, use_fractions=True, time_reversal_convention="spglib")

    centering_symops_xyz = cifblock.get('space_group_symop_magn_centering.xyz')
    if centering_symops_xyz is None:
        centering_symops_alg = cifblock.get('space_group_symop_magn_ssg_centering.algebraic')
        if centering_symops_alg is None:
            centering_symops_xyz = ["x,y,z,+1"]
            cent_symops = _xyzt_symops_to_matrix(
                centering_symops_xyz, use_fractions=True, time_reversal_convention="spglib"
            )
        else:
            centering_symops_xyz = centering_symops_alg
            cent_symops = _alg_symops_to_matrix(
                centering_symops_alg, use_fractions=True, time_reversal_convention="spglib"
            )
    else:
        cent_symops = _xyzt_symops_to_matrix(
            centering_symops_xyz, use_fractions=True, time_reversal_convention="spglib"
        )

    _compose_ops_with_centerings(base_symops, cent_symops)

    bns_nbr = cifblock.get('space_group_magn.number_bns')
    bns_name = cifblock.get('space_group_magn.name_bns')
    parent_name_hm = cifblock.get('parent_space_group.name_h-m_alt')
    space_group_name_hm = re.sub("  +", " ", parent_name_hm.strip()) if parent_name_hm is not None else None
    space_group_nbr = cifblock.get('parent_space_group.it_number')
    icsd = cifblock.get('database_code_ICSD')
    doi = cifblock.get('citation_doi')

    result = {
        'cell_parameters_exact': _cell_parameter_tokens(cifblock),
        'format': 'mcif',
        'positions_exact': asu.positions_exact,
        'occupancies': asu.occupancies,
        'occupancies_exact': asu.occupancies_exact,
        'occupancy_precisions': asu.occupancy_precisions,
        'symbols': asu.symbols,
        'symops_xyz': tuple(raw_symops),
        'incomm': incomm,
        'space_group_nbr': space_group_nbr,
        'space_group_name_hm': space_group_name_hm,
        'icsd': icsd,
        'doi': doi,
        'moment_basis': moment_basis,
        'magmoms_exact': magmoms_exact,
        'centerings_xyz': tuple(centering_symops_xyz),
        'bns_nbr': bns_nbr,
        'bns_name': bns_name,
        'equivalent_atoms': asu.equivalent_atoms,
        'coordinate_precision': asu.coordinate_precision,
        'basis_precision': asu.basis_precision,
        'magmom_precision': magres,
        'labels': asu.labels,
    }

    return result


def mag_asus_from_mcif_file(
    source: str | os.PathLike[str] | Iterable[str], *, error_on_nonmag: bool = False
) -> list[dict[str, Any]]:
    cifblocks, _header = read_cif(source, allow_cif2=True)

    outputs = []
    for name, cifblock in cifblocks:
        outputs += [cifblock_to_mag_asu(cifblock, error_on_nonmag=error_on_nonmag)]
    return outputs


def read_mcif_asus(source: str | os.PathLike[str] | Iterable[str]) -> dict[str, Any]:
    """Read an mcif into the neutral payload used by the ``.mcif`` loader."""
    cifblocks, header = read_cif(source, allow_cif2=True)
    blocks = []
    unparsed = []
    for name, cifblock in cifblocks:
        if 'atom_site_label' not in cifblock:
            continue
        try:
            blocks.append(cifblock_to_mag_asu(cifblock))
        except Exception as error:
            unparsed.append({'block': name, 'reason': f'{type(error).__name__}: {error}'})
    return {'format': 'mcif', 'blocks': blocks, 'unparsed': unparsed, 'header': header}


def single_mag_asu_from_mcif_file(
    source: str | os.PathLike[str] | Iterable[str], *, error_on_nonmag: bool = False
) -> dict[str, Any]:
    cifblocks, _header = read_cif(source, allow_cif2=True)

    # Get the first cifblock with atomic sites
    for name, cifblock in cifblocks:
        if 'atom_site_label' in cifblock:
            break
    else:
        raise ValueError("No structural block found in mcif")

    return cifblock_to_mag_asu(cifblock, error_on_nonmag=error_on_nonmag)
