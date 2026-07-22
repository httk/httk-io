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

import itertools
import math
import re
from fractions import Fraction
from typing import Any

from .cif_parser import (
    parse_asu_cell,
    parse_cif_float,
    parse_linear_expr,
    parse_structural_modulation,
)
from .cif_reader import read_cif


def extract_parent_q_basis(cifblock):
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
        basis.append(tuple(parse_cif_float(v) for v in row))
    return basis


def extract_fourier_coeffs(cifblock, max_q_guess=12):
    """
    Return ``(coeff_rows, m)`` where ``coeff_rows`` is a list of coefficient tuples
    ``(c1, c2, ..., cm)`` and ``m`` is the number of q-vectors detected (``>= 0``).

    All present ``q{i}_coeff`` columns are found, zipped row-wise with missing entries
    filled by zeros, and duplicate coefficient tuples are removed.
    """
    # discover which q*_coeff columns exist
    present_cols = []
    for i in range(1, max_q_guess + 1):
        key = f'atom_site_Fourier_wave_vector.q{i}_coeff'
        col = cifblock.get(key)
        if col is not None:
            present_cols.append((i, key, col))

    if not present_cols:
        return [], 0

    # All present columns must have same length (number of rows). We’ll be permissive: pad shorter ones with zeros.
    max_len = max(len(col) for (_, _, col) in present_cols)
    m = max(i for (i, _, _) in present_cols)

    # Build a dense matrix of size (max_len x m), filling missing cols with zeros
    # 1-based indexing externally; 0-based in list
    columns = [None] * m
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
    rows = list(zip(*columns))

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
        basis = list(basis) + [(0.0, 0.0, 0.0)] * (m - len(basis))
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
    rx, tx = parse_linear_expr(px, use_fractions=use_fractions)
    ry, ty = parse_linear_expr(py, use_fractions=use_fractions)
    rz, tz = parse_linear_expr(pz, use_fractions=use_fractions)
    if time_reversal_convention == "mcif":
        ts = int(ts)
    elif time_reversal_convention == "spglib":
        ts = int((1 - int(ts)) / 2)
    else:
        raise Exception("Unrecognized time reversal convention.")
    return (rx, ry, rz), (tx, ty, tz), ts


def xyzt_symops_to_matrix(symops_xyz, use_fractions=False, time_reversal_convention="mcif"):
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
                raise Exception("Centering symop that includes rotation is invalid")
            t_new = (t[0] + c[0], t[1] + c[1], t[2] + c[2])
            time_new = (time_flag + time_c) % 2  # time_flag * time_c for -1/+1 convention
            composed.append((R, t_new, time_new))
    return composed


def _parse_moments(block, *, k_sigma=2.0, equalize=True, resolution=True) -> tuple[Any, ...] | None:
    """
    Extract magnetic moments from a mcif block.

    Parameters
    ----------
    block: mcif block
    equalize : bool
        If True, perform symmetry-based equalization.
    resolution : bool
        If True, return an additional grid_dens based on numeric resolution.

    Returns
    -------
    If resolution=False:
        moments, labels, spin_basis
    If resolution=True:
        moments, labels, spin_basis, grid_dens
    """

    def _get(name):
        v = block.get(name)
        return [] if v is None else list(v)

    def _len_ok(xs, ys, zs, n):
        return len(xs) == len(ys) == len(zs) == n and n > 0

    labels = _get('atom_site_moment.label')
    n = len(labels)
    if n == 0:
        if resolution:
            return [], [], "crystal", None
        return [], [], "crystal"

    # Try crystal basis first
    xs = _get('atom_site_moment.crystalaxis_x')
    ys = _get('atom_site_moment.crystalaxis_y')
    zs = _get('atom_site_moment.crystalaxis_z')
    spin_basis = "crystal"

    if not _len_ok(xs, ys, zs, n):
        xs = _get('atom_site_moment.Cartn_x')
        ys = _get('atom_site_moment.Cartn_y')
        zs = _get('atom_site_moment.Cartn_z')
        spin_basis = "cartesian"
        if not _len_ok(xs, ys, zs, n):
            return None if not resolution else (None, None, None, None)

    # Labels must be unique
    if n != len(set(labels)):
        raise ValueError("Non-equivalent sites share the same moment label in CIF data.")

    # Read optional equalization metadata only if needed
    if equalize:
        forms = _get('atom_site_moment.symmform')
        mags = _get('atom_site_moment.magnitude')
    else:
        forms = []
        mags = []

    moments = []
    component_resolutions = []

    for i in range(n):
        if resolution:
            mx, mx_meta = parse_cif_float(xs[i], meta=True)
            my, my_meta = parse_cif_float(ys[i], meta=True)
            mz, mz_meta = parse_cif_float(zs[i], meta=True)

            component_resolutions.extend([mx_meta['resolution'], my_meta['resolution'], mz_meta['resolution']])
        else:
            mx = parse_cif_float(xs[i], meta=False)
            my = parse_cif_float(ys[i], meta=False)
            mz = parse_cif_float(zs[i], meta=False)

        # Equalization
        if equalize and i < len(forms):
            form = forms[i].replace(' ', '').lower() if forms[i] else None

            if i < len(mags) and mags[i] not in (None, '?', '.', ''):
                m_val, m_meta = parse_cif_float(mags[i], meta=True)
                m_esd = m_meta.get('esd', None)
            else:
                m_esd = None

            if (
                form in ('mx,mx,mx', 'my,my,my', 'mz,mz,mz')
                and m_esd is not None
                and m_esd > 0.0
                and mx is not None
                and my is not None
                and mz is not None
            ):
                sigma_comp = m_esd / (3.0**0.5)
                mean_comp = (mx + my + mz) / 3.0
                if max(abs(mx - mean_comp), abs(my - mean_comp), abs(mz - mean_comp)) <= k_sigma * sigma_comp:
                    mx = my = mz = mean_comp

        moments.append((mx, my, mz))

    # Fast path: no grid resolution requested
    if not resolution:
        return moments, labels, spin_basis

    # 1. Data resolution = largest implied resolution from written precision
    if component_resolutions:
        data_resolution = max(component_resolutions)
    else:
        data_resolution = 0.0

    # 2. Separation resolution (non-periodic for moments)
    if n > 1:
        diffs = []
        for a, b in itertools.combinations(range(n), 2):
            mx1, my1, mz1 = moments[a]
            mx2, my2, mz2 = moments[b]

            if mx1 is None or my1 is None or mz1 is None or mx2 is None or my2 is None or mz2 is None:
                continue

            for d in (abs(mx1 - mx2), abs(my1 - my2), abs(mz1 - mz2)):
                # Avoid floating noise; only treat real differences as separations
                if d > data_resolution:
                    diffs.append(d)

        if diffs:
            separation_resolution = min(diffs) / 2.0
        else:
            separation_resolution = float('inf')
    else:
        separation_resolution = float('inf')

    if separation_resolution == float('inf'):
        mag_res = data_resolution
    elif data_resolution == 0.0:
        mag_res = separation_resolution
    else:
        mag_res = min(data_resolution, separation_resolution)

    return moments, labels, spin_basis, mag_res


def _parse_linear_expr_algebraic(expr, allowed_vars=('x1', 'x2', 'x3'), use_fractions=False):
    """
    Parse a single algebraic coordinate expression from a superspace op.

    Examples: 'x1-x2', '-x3+1/2', '2x1-x2', 'x1+1/3', 'x1-2x2+3/4'
    Returns (row, const) where:
      row is a list of integer coefficients (may be outside {-1,0,1})
      const is float or Fraction
    Raises if the expression involves variables not in allowed_vars,
    or if a variable has a non-integer coefficient (e.g. 1/2 x1).
    """
    s = expr.replace(" ", "")
    if not s:
        raise ValueError("Empty expression")
    if s[0] not in "+-":
        s = "+" + s

    # ([sign]) ( [optional number] x<digits>  |  standalone number )
    # - number can be: integer, decimal, or fraction a/b
    token_re = (
        r'([+-])'
        r'(?:(?:(?:(\d+(?:/\d+)?|\d*\.\d+)?)'  # optional coefficient before variable
        r'(x\d+))|((?:\d+/\d+)|(?:\d+(?:\.\d+)?)))'  # or standalone number
    )

    coeffs = {v: 0 for v in allowed_vars}
    const = Fraction(0) if use_fractions else 0.0

    pos = 0
    for m in re.finditer(token_re, s):
        if m.start() != pos:
            raise ValueError(f"Unparsed tail in '{expr}' near '{s[pos:]}'")
        pos = m.end()

        sign, coef_str, var, num = m.groups()
        sgn = 1 if sign == '+' else -1

        if var is not None:
            if var not in allowed_vars:
                raise ValueError(f"Expression '{expr}' references {var}, not in allowed {allowed_vars}.")
            # default coefficient is 1 if omitted
            if coef_str in (None, ""):
                coef_val = 1
            else:
                # Coefficients on variables must be integers for symmetry ops
                # Parse via Fraction to catch "2", "2.0", "3/1" etc.
                f = Fraction(coef_str) if '/' in coef_str or '.' in coef_str else Fraction(int(coef_str))
                if f.denominator != 1:
                    raise ValueError(f"Non-integer coefficient {coef_str} on {var} in '{expr}'")
                coef_val = int(f.numerator)
            coeffs[var] += sgn * coef_val
        else:
            # standalone numeric translation
            if use_fractions:
                val = Fraction(num) if '/' in num or '.' in num else Fraction(int(num))
            else:
                val = float(Fraction(num)) if '/' in num else float(num)
            const += sgn * val

    if pos != len(s):
        raise ValueError(f"Unparsed tail in '{expr}' near '{s[pos:]}'")

    const_out = const if use_fractions else float(const)
    row = tuple(int(coeffs[v]) for v in allowed_vars)
    return row, const_out


def parse_alg_op(op, use_fractions=False, time_reversal_convention="mcif"):
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


def alg_symops_to_matrix(symops_alg, use_fractions=False, time_reversal_convention="mcif"):
    return [parse_alg_op(s, use_fractions, time_reversal_convention=time_reversal_convention) for s in symops_alg]


def crystal_to_cartesian(moments_cryst, basis):
    """
    moments_cryst : list of N vectors [mx, my, mz] in crystal coordinates
    basis         : 3x3 list with rows a, b, c in Cartesian
    returns       : list of N vectors in Cartesian coordinates
    """

    def normalize(v):
        n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
        return [v[0] / n, v[1] / n, v[2] / n]

    # --- Extract basis vectors ---

    a = basis[0]
    b = basis[1]
    c = basis[2]

    # --- Unit crystallographic axes ---

    ah = normalize(a)
    bh = normalize(b)
    ch = normalize(c)

    # U matrix columns are ah, bh, ch
    # We only need U^T for multiplication:
    # U^T rows = ah, bh, ch

    result = []

    for m in moments_cryst:
        mx, my, mz = m

        # Row vector multiply by U^T
        x = mx * ah[0] + my * bh[0] + mz * ch[0]
        y = mx * ah[1] + my * bh[1] + mz * ch[1]
        z = mx * ah[2] + my * bh[2] + mz * ch[2]

        result.append([x, y, z])

    return result


def _parse_mag_asu_cell(cifblock, *, moment_equalization=True):
    basis, positions, res, symbols, labels, equivalent_atoms = parse_asu_cell(cifblock)
    moments_result = _parse_moments(cifblock, equalize=moment_equalization, resolution=True)
    assert moments_result is not None  # resolution=True always yields a 4-tuple
    cif_moments, momlabels, spin_basis, magres = moments_result

    if spin_basis == "crystal":
        cif_moments = crystal_to_cartesian(cif_moments, basis)
        spin_basis = "cartesian"

    moments_map = {label: (mom[0], mom[1], mom[2]) for label, mom in zip(momlabels, cif_moments)}
    magmoms = [moments_map[i] if i in moments_map else (0.0, 0.0, 0.0) for i in labels]

    return basis, positions, res, magmoms, spin_basis, magres, symbols, labels, equivalent_atoms


def _get_magnetic_fourier_info(cifblock):
    has = False
    atoms = set()

    labels = cifblock.get('_atom_site_moment_Fourier.atom_site_label')
    if labels:
        has = True
        atoms.update(labels)

    return has, sorted(atoms)


def _parse_modulation(cifblock):
    structural_q, mod_dim, has_struct_mod, struct_mod_atoms = parse_structural_modulation(cifblock)

    # Magnetic q
    magnetic_q = None

    # (A) magnetic superspace -> uses same q as structural superspace
    if '_space_group.magn_ssg_name' in cifblock and structural_q:
        magnetic_q = structural_q

    # (B) commensurate magnetic propagation vector
    elif cifblock.get('parent_propagation_vector.kxkykz'):
        rows = cifblock['parent_propagation_vector.kxkykz']
        magnetic_q = [[parse_cif_float(v) for v in row] for row in rows]

    # (C) Fourier-defined magnetic propagation vector
    else:
        fourier = extract_fourier(cifblock)
        if fourier:
            magnetic_q = _collect_k_from_fourier(fourier)

    has_mag_mod, mag_mod_atoms = _get_magnetic_fourier_info(cifblock)

    return structural_q, magnetic_q, mod_dim, has_struct_mod, has_mag_mod, struct_mod_atoms, mag_mod_atoms


def is_rational_component(x, max_den=12, tol=1e-6):
    from fractions import Fraction

    fx = Fraction(x).limit_denominator(max_den)
    return abs(float(fx) - x) < tol


def cifblock_to_mag_asu(cifblock, *, error_on_nonmag=False):

    basis, positions, res, magmoms, spin_basis, magres, symbols, labels, equivalent_atoms = _parse_mag_asu_cell(
        cifblock
    )
    structural_q, magnetic_q, mod_dim, has_struct_mod, has_mag_mod, struct_mod_atoms, mag_mod_atoms = _parse_modulation(
        cifblock
    )

    if error_on_nonmag and magmoms is None:
        raise Exception("Could not extract magnetic moments from mcif file")

    # Determine if magnetic q is incommensurate
    mq_is_incomm = False
    if magnetic_q:
        for q in magnetic_q:
            if any(not is_rational_component(x) for x in q):
                mq_is_incomm = True
                break

    # Build incommensurate descriptor only when really needed
    incomm = None
    if mod_dim > 0 or mq_is_incomm:
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
            raise Exception("No symmetry operations in mcif")
        base_symops = alg_symops_to_matrix(base_symops_alg, use_fractions=True, time_reversal_convention="spglib")
    else:
        base_symops = xyzt_symops_to_matrix(base_symops_xyz, use_fractions=True, time_reversal_convention="spglib")

    centering_symops_xyz = cifblock.get('space_group_symop_magn_centering.xyz')
    if centering_symops_xyz is None:
        centering_symops_alg = cifblock.get('space_group_symop_magn_ssg_centering.algebraic')
        if centering_symops_alg is None:
            centering_symops_xyz = ["x,y,z,+1"]
            cent_symops = xyzt_symops_to_matrix(
                centering_symops_xyz, use_fractions=True, time_reversal_convention="spglib"
            )
        else:
            cent_symops = alg_symops_to_matrix(
                centering_symops_alg, use_fractions=True, time_reversal_convention="spglib"
            )
    else:
        cent_symops = xyzt_symops_to_matrix(centering_symops_xyz, use_fractions=True, time_reversal_convention="spglib")

    symops = _compose_ops_with_centerings(base_symops, cent_symops)

    bns_nbr = cifblock.get('space_group_magn.number_bns')
    bns_name = cifblock.get('space_group_magn.name_bns')
    space_group_name_hm = re.sub("  +", " ", cifblock.get('parent_space_group.name_h-m_alt').strip())
    space_group_nbr = cifblock.get('parent_space_group.it_number')
    icsd = cifblock.get('database_code_ICSD')
    doi = cifblock.get('citation_doi')

    result = {
        'basis': basis,
        'positions': positions,
        'symbols': symbols,
        'symops': symops,
        'incomm': incomm,
        'space_group_nbr': space_group_nbr,
        'space_group_name_hm': space_group_name_hm,
        'icsd': icsd,
        'doi': doi,
        'spin_basis': spin_basis,
        'magmoms': magmoms,
        'bns_nbr': bns_nbr,
        'bns_name': bns_name,
        'equivalent_atoms': equivalent_atoms,
        'resolution': res,
        'magmom_resolution': magres,
        'labels': labels,
    }

    return result


def mag_asus_from_mcif_file(fs, *, error_on_nonmag=False):
    cifblocks, header = read_cif(fs, allow_cif2=True)

    outputs = []
    for name, cifblock in cifblocks:
        outputs += [cifblock_to_mag_asu(cifblock, error_on_nonmag=error_on_nonmag)]
    return outputs


def single_mag_asu_from_mcif_file(fs, *, error_on_nonmag=False):
    cifblocks, header = read_cif(fs, allow_cif2=True)

    # Get the first cifblock with atomic sites
    for name, cifblock in cifblocks:
        if 'atom_site_label' in cifblock:
            break
    else:
        raise Exception("No structural block found in CIF.")

    return cifblock_to_mag_asu(cifblock, error_on_nonmag=error_on_nonmag)
