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
from decimal import Decimal
from fractions import Fraction
from typing import Any, Literal, TypedDict, overload

from .cif_reader import read_cif


class CifMeta(TypedDict):
    """Metadata returned by :func:`parse_cif_float` when ``meta=True``."""

    esd: float | None
    resolution: float


# Regexp close to https://www.iucr.org/__data/iucr/cifdic_html/2/cif_mm.dic/Dtypecodes.html
# matches:  1.234(5), -12.3(12), 3(1)E2, 1.0e-3, +4.2, etc.
_CIF_NUM_RE = re.compile(
    r'^(?P<sign>-)?'  # optional leading minus
    r'(?P<mant>(?:\d+\.?|\d*\.\d+))(\((?P<esd>\d+)\))?'  # mantissa + optional (uncertainty)
    r'(?:[eE](?P<exp>[+-]?\d+))?$'  # optional exponent
)


def _literal_resolution(txt: str) -> float:
    """
    Resolution implied by the literal format of txt (no ESD logic).

    - '1/3'      -> 0.0  (treated as exact)
    - '0.123'    -> 1e-3
    - '5.'       -> 1.0
    - '.25'      -> 1e-2
    - '10'       -> 1.0
    """
    if txt is None:
        return 0.0
    s = txt.strip()
    if not s:
        return 0.0

    # Strip leading sign
    if s[0] in "+-":
        s = s[1:]

    # Fractions -> exact
    if "/" in s:
        try:
            Fraction(s)  # Only to validate it
            return 0.0
        except Exception:
            return 0.0

    # No exponent here if we use mant_str; but be robust
    s = s.lower()
    if "e" in s:
        s, _ = s.split("e", 1)

    if "." in s:
        before, after = s.split(".", 1)
        # malformed: just ignore non-digits
        digits = "".join(ch for ch in after if ch.isdigit())
        if not digits:  # '5.' or '.' or '.e3'
            return 1.0
        return 10.0 ** (-len(digits))
    else:
        # Integer literal
        return 1.0


@overload
def parse_cif_float(token: str, *, meta: Literal[False] = ..., pragmatic: bool = ...) -> float | None: ...


@overload
def parse_cif_float(token: str, *, meta: Literal[True], pragmatic: bool = ...) -> tuple[float | None, CifMeta]: ...


def parse_cif_float(
    token: str, *, meta: bool = False, pragmatic: bool = False
) -> float | None | tuple[float | None, CifMeta]:
    """
    Parse a CIF numeric field.

    If meta=False:
        return float_value or None.

    If meta=True:
        return (float_value_or_None, {'esd': esd_or_None, 'resolution': float_resolution}).
    """
    if token is None:
        raise Exception("parse_cif_float parsing None")

    t = token.strip()
    if t == '?':
        if meta:
            return None, {'esd': None, 'resolution': 0.0}
        return None

    if t in ('.', ''):
        raise Exception("Missing cif value cannot be conveted to float")

    # Replace unicode minus
    if any(ch in t for ch in ("\u2212", "\u2013", "\u2014")):
        if pragmatic:
            t = t.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
        else:
            raise Exception("Cif contains non-ascii minus sign: " + str(t))

    m = _CIF_NUM_RE.match(t)

    if not m:
        # fractions are allowed here too
        try:
            if "/" in t:
                val = float(Fraction(t))
            else:
                val = float(t)
        except Exception:
            # last resort: grab first float-looking chunk
            val = float(re.split(r'([0-9]*(\.[0-9]+)?)', t)[1])

        res = _literal_resolution(t)
        if meta:
            return val, {'esd': None, 'resolution': res}
        return val

    # ---- Normal CIF number ----
    sign = -1 if m.group('sign') == '-' else 1
    mant_str = m.group('mant')
    mant = Decimal(mant_str)
    exp = int(m.group('exp') or '0')
    val = float(sign * mant * (Decimal(10) ** exp))

    # resolution from mantissa literal
    res = _literal_resolution(mant_str) * (10.0**exp)

    esd_str = m.group('esd')
    if not meta:
        # Ignore esd; user didn't ask for meta info
        return val

    if esd_str is not None:
        # Classic CIF esd logic
        if '.' in mant_str:
            dec_places = len(mant_str.split('.', 1)[1])
        else:
            dec_places = 0
        esd_abs = Decimal(int(esd_str)) * (Decimal(10) ** (exp - dec_places))
        esd_val = float(esd_abs)
    else:
        esd_val = None

    return val, {'esd': esd_val, 'resolution': res}


def parse_cif_int(token: str, *, strict: bool = True, allow_round: bool = False) -> int:
    """
    Convert a CIF numeric token (e.g., '123(4)', '3E2', '1.0E3') to an int using the central value.
    - strict=True: require the value to be exactly integral; otherwise raise ValueError.
    - allow_round=True (only if strict=False): round half-even to the nearest int.
    """
    t = token.strip()
    if t in ('.', '?', ''):
        raise ValueError("Missing CIF value cannot be converted to int")

    m = _CIF_NUM_RE.match(t)
    if not m:
        # Fall back for plain integers without (esd)/exponent; will raise if not int-like
        val = Decimal(t)
    else:
        sign = -1 if m.group('sign') == '-' else 1
        mant = Decimal(m.group('mant'))
        exp = int(m.group('exp') or '0')
        val = sign * mant * (Decimal(10) ** exp)

    # Decide how to coerce
    if strict:
        # exactly integral?
        if val == val.to_integral_value():  # no fractional part
            return int(val)
        raise ValueError(f"Non-integer numeric cannot be coerced strictly: {token!r}")
    else:
        # Non-strict: either require integral or allow rounding
        if val == val.to_integral_value():
            return int(val)
        if allow_round:
            return int(val.to_integral_value())  # banker's rounding (half-even)
        raise ValueError(f"Non-integer numeric (set allow_round=True to round): {token!r}")


def parse_linear_expr(expr, use_fractions=False):
    """
    expr: e.g. 'x-y', '-z+1/2', 'x', 'y', 'z-1', 'x-2y', '3x+1/2'
    Returns (row, const) where row is [ax, ay, az] (integers or Fractions),
    and const is a float or Fraction depending on use_fractions.
    """
    s = expr.replace(" ", "")
    if not s:
        raise ValueError("Empty expression")
    if s[0] not in "+-":
        s = "+" + s

    coeffs: dict[str, Fraction | float | int] = {'x': 0, 'y': 0, 'z': 0}
    const = Fraction(0) if use_fractions else 0.0

    # ([sign]) ( [optional number] [var]  |  number )
    token_re = r'([+-])(?:(?:(?:(\d+(?:/\d+)?|\d*\.\d+)?)' r'(x|y|z))|((?:\d+/\d+)|(?:\d+(?:\.\d+)?)))'

    pos = 0
    for m in re.finditer(token_re, s):
        if m.start() != pos:
            # There are leftover characters -> invalid tokenization
            raise ValueError(f"Unparsed tail in '{expr}' near '{s[pos:]}'")
        pos = m.end()

        sign, coef_str, var, num = m.groups()
        sgn = 1 if sign == '+' else -1

        if var is not None:
            # variable term ± (coef or 1) * var
            if coef_str in (None, ""):
                coef_val = 1
            else:
                coef_val = Fraction(coef_str) if '/' in coef_str else float(coef_str)
            if use_fractions and not isinstance(coef_val, Fraction):
                coef_val = Fraction(coef_val)
            coeffs[var] += sgn * coef_val
        else:
            # standalone numeric translation
            if use_fractions:
                val = Fraction(num) if '/' in num else Fraction(num)
            else:
                val = float(Fraction(num)) if '/' in num else float(num)
            const += sgn * val

    if pos != len(s):
        raise ValueError(f"Unparsed tail in '{expr}' near '{s[pos:]}'")

    # Ensure output constant has requested type
    const_out = const if use_fractions else float(const)

    # Order as (x,y,z)
    return (coeffs['x'], coeffs['y'], coeffs['z']), const_out


def parse_xyz_op(op, use_fractions=False):
    """
    op: e.g. 'x-y,x,-z+1/2'
    Returns (R, t) where R is 3x3 (list of rows), t is length-3 float list
    """
    parts = [p.strip() for p in op.split(",")]
    if len(parts) != 3:
        raise ValueError(f"Unexpected op format: {op}")
    px, py, pz = parts
    rx, tx = parse_linear_expr(px, use_fractions=use_fractions)
    ry, ty = parse_linear_expr(py, use_fractions=use_fractions)
    rz, tz = parse_linear_expr(pz, use_fractions=use_fractions)
    return (rx, ry, rz), (tx, ty, tz)


def xyz_symops_to_matrix(symops_xyz, use_fractions=False):
    return [parse_xyz_op(s, use_fractions) for s in symops_xyz]


def cif_exact_token(token: str) -> str | None:
    """The numeric part of a CIF value, as text, with any uncertainty estimate removed.

    ``"0.3333(5)"`` becomes ``"0.3333"``; ``"?"`` and ``"."`` become ``None``. The point of
    keeping the text rather than a float is fidelity: a consumer that wants an exact value
    can read ``0.3333`` as the rational 3333/10000, which is what the file says, whereas
    ``float("0.3333")`` is a binary approximation whose exact rational value is
    6004199023210345/18014398509481984 and says something the file did not.
    """
    text = token.strip().strip("'\"")
    if text in ("?", ".", ""):
        return None
    if "(" in text:
        text = text[: text.index("(")]
    return text.strip() or None


def _parse_atoms(block, resolution=True) -> tuple[Any, ...]:
    """
    Returns:
      if resolution == False:
         (symbols, labels, positions, exact_positions, occupancies)
      if resolution == True:
         (symbols, labels, positions, exact_positions, occupancies, res)

    positions: list of (x, y, z) floats
    exact_positions: list of (x, y, z) numeric strings, uncertainties stripped, so a
      consumer can embed the coordinate the file actually wrote rather than its binary
      float approximation
    occupancies:
      - list of floats (same length as symbols) if '_atom_site_occupancy' exists
      - None otherwise
    """
    syms = block.get('atom_site_type_symbol')
    lbs = block.get('atom_site_label')
    xs = block.get('atom_site_fract_x')
    ys = block.get('atom_site_fract_y')
    zs = block.get('atom_site_fract_z')

    n = len(xs)
    assert len(ys) == len(zs) == len(lbs) == len(syms) == n

    # Optional occupancy column
    occ_col = block.get('atom_site_occupancy')
    if occ_col is not None:
        occs = []
        for t in occ_col:
            v = parse_cif_float(t, meta=False)
            occs.append(v)
    else:
        occs = None

    symbols = [s.strip() for s in syms]
    labels = [lab.strip() for lab in lbs]

    exact_positions = [
        (cif_exact_token(xi), cif_exact_token(yi), cif_exact_token(zi)) for xi, yi, zi in zip(xs, ys, zs)
    ]

    # Fast path: no resolution / grid requested
    if not resolution:
        positions = [
            (
                parse_cif_float(xi, meta=False),
                parse_cif_float(yi, meta=False),
                parse_cif_float(zi, meta=False),
            )
            for xi, yi, zi in zip(xs, ys, zs)
        ]
        return symbols, labels, positions, exact_positions, occs

    # Full path with resolution
    positions = []
    coord_resolutions = []

    for xi, yi, zi in zip(xs, ys, zs):
        vx, mx = parse_cif_float(xi, meta=True)
        vy, my = parse_cif_float(yi, meta=True)
        vz, mz = parse_cif_float(zi, meta=True)

        positions.append((vx, vy, vz))
        coord_resolutions.extend(
            [
                mx.get('resolution', 0.0),
                my.get('resolution', 0.0),
                mz.get('resolution', 0.0),
            ]
        )

    # 1) Data-implied resolution: take the COARSEST (largest) non-zero resolution.
    finite_res = [r for r in coord_resolutions if r is not None]
    if finite_res:
        data_resolution = max(finite_res)
    else:
        data_resolution = 0.0  # no constraint from data formatting

    # 2) Separation resolution from coordinates (fractional, so use periodic deltas)
    if n > 1:

        def periodic_delta(a, b):
            d = abs(a - b)
            return min(d, 1.0 - d)  # fractional periodicity

        eps = data_resolution / 10 if data_resolution > 0.0 else 1e-12

        # compute deltas for each axis
        deltas_x = [periodic_delta(positions[i][0], positions[j][0]) for i, j in itertools.combinations(range(n), 2)]
        deltas_y = [periodic_delta(positions[i][1], positions[j][1]) for i, j in itertools.combinations(range(n), 2)]
        deltas_z = [periodic_delta(positions[i][2], positions[j][2]) for i, j in itertools.combinations(range(n), 2)]

        # filter out zero-ish separations using adaptive epsilon
        pos_x = [d for d in deltas_x if d > eps]
        pos_y = [d for d in deltas_y if d > eps]
        pos_z = [d for d in deltas_z if d > eps]

        sep_x = min(pos_x) if pos_x else float('inf')
        sep_y = min(pos_y) if pos_y else float('inf')
        sep_z = min(pos_z) if pos_z else float('inf')

        min_sep = min(sep_x, sep_y, sep_z)

        if math.isinf(min_sep):
            separation_resolution = float('inf')
        else:
            separation_resolution = min_sep / 2.0

    else:
        separation_resolution = float('inf')

    # 3) Final res
    if separation_resolution == float('inf'):
        res = data_resolution
    elif data_resolution == 0.0:
        res = separation_resolution
    else:
        res = min(data_resolution, separation_resolution)

    # grid should sit in [0, 1], but fractional coords guarantee that anyway.
    return symbols, labels, positions, exact_positions, occs, res


def _parse_uc(block):
    """The six cell parameters, or a message naming the tags the block is missing.

    Reporting the missing tags matters because an incomplete cell is the most common
    reason a data block turns out not to be a structure, and "parse_cif_float parsing
    None" tells a reader nothing about which tag to go and look for.
    """
    tags = (
        'cell_length_a',
        'cell_length_b',
        'cell_length_c',
        'cell_angle_alpha',
        'cell_angle_beta',
        'cell_angle_gamma',
    )
    missing = [tag for tag in tags if block.get(tag) is None]
    if missing:
        raise ValueError(f"CIF block has no unit cell: missing {', '.join('_' + tag for tag in missing)}")
    return tuple(parse_cif_float(block.get(tag)) for tag in tags)


def _basis_from_lengths_angles(a, b, c, alpha, beta, gamma):
    """
    Conventional 3x3 lattice (rows are a,b,c in Cartesian Å) from a,b,c (Å) and angles (deg).
    """

    def _deg2rad(d):
        return d * math.pi / 180.0

    alpha, beta, gamma = map(_deg2rad, (alpha, beta, gamma))
    ca, cb, cg = math.cos(alpha), math.cos(beta), math.cos(gamma)
    sg = math.sin(gamma)

    ax, ay, az = a, 0.0, 0.0
    bx, by, bz = b * cg, b * sg, 0.0
    # cz via the standard formula for triclinic cells
    cx = c * cb
    cy = c * (ca - cb * cg) / (sg if abs(sg) > 1e-12 else 1.0)
    cz_sq = c**2 - cx**2 - cy**2
    cz = math.sqrt(max(cz_sq, 0.0))
    return [[ax, ay, az], [bx, by, bz], [cx, cy, cz]]


def parse_asu_cell(cifblock):
    a, b, c, alpha, beta, gamma = _parse_uc(cifblock)
    basis = _basis_from_lengths_angles(a, b, c, alpha, beta, gamma)
    symbols, labels, positions, exact_positions, occs, res = _parse_atoms(cifblock, resolution=True)

    # figure out equivalent atoms based on labels
    labels_map = {}
    equivalent_atoms = []
    next_id = 1
    for lab in labels:
        if lab not in labels_map:
            labels_map[lab] = next_id
            next_id += 1
        equivalent_atoms.append(labels_map[lab])

    return basis, positions, exact_positions, occs, res, symbols, labels, equivalent_atoms


def parse_structural_modulation(cifblock):
    """
    Extract structural superspace modulation information from a standard CIF.

    Returns a tuple ``(structural_q, mod_dim, has_struct_mod, struct_mod_atoms)`` where
    ``structural_q`` is a list of q-vectors or ``None``, ``mod_dim`` is the modulation
    dimension (0 if absent), ``has_struct_mod`` is a bool, and ``struct_mod_atoms`` is a
    sorted list of atom-site labels.
    """
    # modulation dimension (0 if absent)
    mod_dim = int(cifblock.get('cell_modulation_dimension', 0))

    # structural_q from cell_wave_vector (only if mod_dim > 0)
    structural_q = None
    qx = cifblock.get('_cell_wave_vector_x')
    qy = cifblock.get('_cell_wave_vector_y')
    qz = cifblock.get('_cell_wave_vector_z')
    if qx and qy and qz:
        structural_q = [[float(qx[i]), float(qy[i]), float(qz[i])] for i in range(len(qx))]

    # detect structural Fourier modulations
    has_struct_mod = False
    struct_mod_atoms = set()

    labels = cifblock.get('_atom_site_displace_Fourier.atom_site_label')
    if labels:
        has_struct_mod = True
        struct_mod_atoms.update(labels)

    labels = cifblock.get('_atom_site_occupancy_Fourier.atom_site_label')
    if labels:
        has_struct_mod = True
        struct_mod_atoms.update(labels)

    return structural_q, mod_dim, has_struct_mod, sorted(struct_mod_atoms)


def _first_tag(cifblock, *names):
    """The value of the first tag present, or ``None``.

    Tags are matched against the reader's normalized spelling (lower case, no leading
    underscore) and, defensively, against the name as written, so a differently normalizing
    reader would still resolve them.
    """
    for name in names:
        for candidate in (name, name.lower(), '_' + name):
            value = cifblock.get(candidate)
            if value is not None:
                return value
    return None


def cifblock_to_asu(cifblock, *, return_single=False):

    # basic atom-site parsing
    (
        basis,
        positions,
        exact_positions,
        occupancies,
        resolution,
        symbols,
        labels,
        equivalent_atoms,
    ) = parse_asu_cell(cifblock)

    # standard space group symmetry
    symops_xyz = cifblock.get('space_group_symop.operation_xyz')
    if symops_xyz is None:
        # Some readers normalize loop tags without dots, so accept that too.
        symops_xyz = cifblock.get('space_group_symop_operation_xyz')
    if symops_xyz is None:
        # some CIFs use older spelling
        symops_xyz = cifblock.get('symmetry_equiv_pos_as_xyz')

    if symops_xyz is None:
        raise Exception("No symmetry operations in CIF.")

    symops = xyz_symops_to_matrix(symops_xyz, use_fractions=True)

    # structural modulation
    structural_q, mod_dim, has_struct_mod, struct_atoms = parse_structural_modulation(cifblock)

    # Build the incommensurate structure descriptor, or None
    incomm = None
    if mod_dim > 0 or structural_q or has_struct_mod:
        incomm = {
            'structural_q': structural_q,
            'mod_dim': mod_dim,
            'has_structural_modulation': has_struct_mod,
            'structural_modulated_atoms': struct_atoms,
        }

    # The reader normalizes tags to lower case with the leading underscore stripped, so
    # these lookups must be spelled that way; written with the underscore and the original
    # capitalization they silently matched nothing and every CIF came back with no space
    # group at all.
    space_group_name_hm = _first_tag(cifblock, 'space_group_name_h-m_alt', 'symmetry_space_group_name_h-m')
    space_group_name_hall = _first_tag(cifblock, 'space_group_name_hall', 'symmetry_space_group_name_hall')
    space_group_nbr = _first_tag(cifblock, 'space_group_it_number', 'symmetry_space_group_it_number')
    icsd = _first_tag(cifblock, 'database_code_icsd')
    doi = _first_tag(cifblock, 'citation_doi')

    return {
        'format': 'cif',
        'basis': basis,
        'cell_parameters': _parse_uc(cifblock),
        'positions': positions,
        'positions_exact': exact_positions,
        'occupancies': occupancies,
        'symbols': symbols,
        'symops': symops,
        'incomm': incomm,
        'space_group_nbr': space_group_nbr,
        'space_group_name_hm': space_group_name_hm,
        'space_group_name_hall': space_group_name_hall,
        'icsd': icsd,
        'doi': doi,
        'resolution': resolution,
        'equivalent_atoms': equivalent_atoms,
        'labels': labels,
    }


def read_cif_asus(fs):
    """Read a CIF and return its asymmetric units as a neutral, tagged payload.

    This is what ``httk.core.load`` returns for a ``.cif`` file: a mapping with
    ``format`` set to ``"cif"``, ``blocks`` holding one asymmetric-unit mapping per data
    block that describes a structure, and ``header`` the file's leading comment. The tag
    is what lets a consumer dispatch on the file type without knowing which reader
    produced the payload, the same way the POSCAR reader already works.

    Loading never fails on account of a block that is not a structure. CIF is a
    general-purpose format and a file may hold bibliographic entries, powder patterns, or
    an incomplete draft alongside — or instead of — anything crystallographic. Blocks
    without atom sites are simply not structures and are passed over; blocks that have
    atom sites but cannot be interpreted are collected in ``unparsed``, each with the
    reason, so that nothing is dropped silently and the failure surfaces when a structure
    is actually asked for.

    The mapping stays neutral — plain lists, strings, and exact ``Fraction`` symmetry
    operations, with no domain objects — so *httk-io* need not know about
    *httk-atomistic*. Turning it into a structure is
    ``httk.atomistic.asu_structure_from_cif``.
    """
    cifblocks, header = read_cif(fs, allow_cif2=False)

    blocks = []
    unparsed = []
    for name, cifblock in cifblocks:
        if 'atom_site_label' not in cifblock:
            continue
        try:
            blocks.append(cifblock_to_asu(cifblock))
        except Exception as error:
            unparsed.append({'block': name, 'reason': f'{type(error).__name__}: {error}'})

    return {'format': 'cif', 'blocks': blocks, 'unparsed': unparsed, 'header': header}


def asus_from_cif_file(fs):
    cifblocks, header = read_cif(fs, allow_cif2=False)

    outputs = []
    for name, cifblock in cifblocks:
        outputs += [cifblock_to_asu(cifblock)]
    return outputs


def single_asu_from_cif_file(fs):
    cifblocks, header = read_cif(fs, allow_cif2=False)

    # Get the first cifblock with atomic sites
    for name, cifblock in cifblocks:
        if 'atom_site_label' in cifblock:
            break
    else:
        raise Exception("No structural block found in CIF.")

    return cifblock_to_asu(cifblock)
