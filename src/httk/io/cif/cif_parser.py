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

import math
import os
import re
import warnings
from collections.abc import Iterable, Mapping
from decimal import Decimal
from fractions import Fraction
from typing import Any, Literal, NamedTuple, TypedDict, cast, overload

from httk.core import combined_precision, decimal_precision

from .cif_reader import read_cif
from .cif_tags import CIF_TAGS


class CifMeta(TypedDict):
    """What a CIF number claims about its own precision, beyond its value.

    Both are exact rationals, or ``None`` when the token makes no such claim. ``None`` is
    not the same as claiming exactness: a value written ``1/3`` states an exact number and
    reports ``None`` for both, while ``?`` states nothing at all and also reports ``None``.

    ``precision`` is implied by the digits written (``0.3333`` is good to ``1/10000``);
    ``esd`` is the standard uncertainty a file states explicitly, so ``0.3333(7)`` reports
    a precision of ``1/10000`` *and* an esd of ``7/10000``. They measure different things
    and the coarser of the two is the honest claim.

    :param esd: Standard uncertainty explicitly stated by the CIF token.
    :param precision: Precision implied by the digits written in the CIF token.
    """

    esd: Fraction | None
    precision: Fraction | None


# IUPAC element symbols (H..Og); used to infer a site's element from its label when a CIF
# omits the optional _atom_site_type_symbol column.
_ELEMENT_SYMBOLS = frozenset(
    (
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
        'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
        'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
        'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
        'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
        'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
        'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
        'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
        'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
    )
)  # fmt: skip

# The leading alphabetic run of an atom-site label, e.g. "Mg" of "MgM1".
_LABEL_PREFIX_RE = re.compile(r'[A-Za-z]+')

# Regexp close to https://www.iucr.org/__data/iucr/cifdic_html/2/cif_mm.dic/Dtypecodes.html
# matches:  1.234(5), -12.3(12), 3(1)E2, 1.0e-3, +4.2, etc.
_CIF_NUM_RE = re.compile(
    r'^(?P<sign>[+-])?'  # optional leading sign
    r'(?P<mant>(?:\d+\.?|\d*\.\d+))(\((?P<esd>\d+)\))?'  # mantissa + optional (uncertainty)
    r'(?:[eE](?P<exp>[+-]?\d+))?$'  # optional exponent
)


@overload
def parse_cif_float(token: str, *, meta: Literal[False] = ..., pragmatic: bool = ...) -> float | None: ...


@overload
def parse_cif_float(token: str, *, meta: Literal[True], pragmatic: bool = ...) -> tuple[float | None, CifMeta]: ...


def parse_cif_float(
    token: str, *, meta: bool = False, pragmatic: bool = False
) -> float | None | tuple[float | None, CifMeta]:
    """Parse a CIF numeric field.

    If meta=False:
        return float_value or None.

    If meta=True:
        return ``(float_value_or_None, CifMeta)`` — the value plus what the token claims
        about its own precision. See :class:`CifMeta`; both of its fields are exact
        rationals or ``None``.

    :param token: Numeric text to parse, including CIF uncertainty notation when present.
    :param meta: Include the value's claimed precision and standard uncertainty.
    :param pragmatic: Salvage selected non-standard numeric spellings instead of rejecting them.
    :return: The parsed value, optionally paired with its precision metadata.
    :raises ValueError: If the token is missing or cannot be interpreted as a CIF number.
    """
    if token is None:
        raise ValueError("Cannot parse None as a CIF float")

    t = token.strip()
    if t == '?':
        if meta:
            return None, {'esd': None, 'precision': None}
        return None

    if t in ('.', ''):
        raise ValueError("Missing CIF value cannot be converted to float")

    # Replace unicode minus
    if any(ch in t for ch in ("\u2212", "\u2013", "\u2014")):
        if pragmatic:
            t = t.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
        else:
            raise ValueError("CIF contains a non-ASCII minus sign: " + str(t))

    m = _CIF_NUM_RE.match(t)

    if not m:
        # fractions are allowed here too
        try:
            if "/" in t:
                val = float(Fraction(t))
            else:
                val = float(t)
        except (ValueError, ZeroDivisionError) as error:
            if not pragmatic:
                raise ValueError(f"Invalid CIF numeric token: {token!r}") from error
            warnings.warn(f"Salvaging malformed CIF numeric token: {token!r}", RuntimeWarning, stacklevel=2)
            try:
                val = float(re.split(r'([0-9]*(\.[0-9]+)?)', t)[1])
            except (IndexError, ValueError) as salvage_error:
                raise ValueError(f"Invalid CIF numeric token: {token!r}") from salvage_error

        if meta:
            # Anything the CIF number pattern did not match: a fraction, or salvage. A
            # fraction states an exact value, and salvaged junk states nothing reliable,
            # so decimal_precision returns None for both and that is the honest answer.
            return val, {'esd': None, 'precision': decimal_precision(t)}
        return val

    # ---- Normal CIF number ----
    sign = -1 if m.group('sign') == '-' else 1
    mant_str = m.group('mant')
    mant = Decimal(mant_str)
    exp = int(m.group('exp') or '0')
    val = float(sign * mant * (Decimal(10) ** exp))

    # Precision implied by the digits written, scaled by any exponent.
    precision = decimal_precision(mant_str)
    if precision is not None and exp:
        precision = precision * Fraction(10) ** exp

    esd_str = m.group('esd')
    if not meta:
        # Ignore esd; user didn't ask for meta info
        return val

    if esd_str is not None:
        # Classic CIF esd logic: the parenthesised integer applies to the last written
        # digit, so 1.234(5) is 1.234 +/- 0.005. Kept exact rather than rounded to a float.
        if '.' in mant_str:
            dec_places = len(mant_str.split('.', 1)[1])
        else:
            dec_places = 0
        esd_val: Fraction | None = Fraction(int(esd_str)) * Fraction(10) ** (exp - dec_places)
    else:
        esd_val = None

    return val, {'esd': esd_val, 'precision': precision}


def parse_cif_fraction(token: str) -> Fraction | None:
    """Parse a CIF numeric token exactly, preserving finite decimals and fractions.

    :param token: Numeric text whose exact central value should be preserved.
    :return: The exact value, or ``None`` for an unknown CIF value.
    :raises ValueError: If the token does not contain a valid exact numeric value.
    """
    exact = cif_exact_token(token)
    if exact is None:
        return None
    try:
        return Fraction(exact)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError(f"Invalid exact CIF numeric token: {token!r}") from error


def parse_cif_int(token: str, *, strict: bool = True, allow_round: bool = False) -> int:
    """Convert a CIF numeric token to an integer using its central value.

    The accepted forms include ``'123(4)'``, ``'3E2'``, and ``'1.0E3'``.
    - strict=True: require the value to be exactly integral; otherwise raise ValueError.
    - allow_round=True (only if strict=False): round half-even to the nearest int.

    :param token: Numeric text to convert.
    :param strict: Require the central value to be integral.
    :param allow_round: Permit half-even rounding when ``strict`` is false.
    :return: The converted integer.
    :raises ValueError: If the token is missing or its central value cannot be converted under the selected rules.
    :raises decimal.InvalidOperation: If an unmatched token has malformed decimal syntax.
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


def cif_exact_token(token: str) -> str | None:
    """The numeric part of a CIF value, as text, with any uncertainty estimate removed.

    ``"0.3333(5)"`` becomes ``"0.3333"`` and ``"3(1)e-1"`` becomes ``"3e-1"``;
    ``"?"`` and ``"."`` become ``None``. The point of keeping the text rather than a
    float is fidelity: a consumer that wants an exact value can read ``0.3333`` as the
    rational 3333/10000, which is what the file says, whereas ``float("0.3333")`` is a
    binary approximation whose exact rational value is 6004199023210345/18014398509481984
    and says something the file did not.

    :param token: CIF value text whose uncertainty estimate should be removed.
    :return: The central numeric text, or ``None`` for an unknown value.
    """
    text = token.strip().strip("'\"")
    if text in ("?", ".", ""):
        return None
    # The exponent belongs to the central value and follows the parenthesized ESD in CIF
    # syntax, so remove only the ESD instead of truncating everything after its opening
    # parenthesis.
    text = re.sub(r'\(\d+\)(?=(?:[eE][+-]?\d+)?$)', '', text)
    return text.strip() or None


def _symbol_from_label(label: str) -> str | None:
    """Element symbol inferred from an ``_atom_site_label``, or ``None``.

    :param label: An atom-site label such as ``"MgM1"`` or ``"O1"``.
    :return: The inferred IUPAC element symbol, or ``None`` when no leading element is recognised.
    """
    match = _LABEL_PREFIX_RE.match(label.strip())
    if match is None:
        return None
    prefix = match.group(0)
    if len(prefix) >= 2:
        # A two-letter symbol only when the case supports it: an "Xy" prefix reads as one,
        # and an all-caps label ("MGM1") still resolves, but a mixed-case "OSi1" must stay
        # "O" rather than becoming the greedy two-letter "Os".
        two = prefix[:2].capitalize()
        if two in _ELEMENT_SYMBOLS and (prefix[1].islower() or prefix.isupper()):
            return two
    one = prefix[:1].upper()
    if one in _ELEMENT_SYMBOLS:
        return one
    return None


def _parse_atoms(block: Mapping[str, Any]) -> tuple[Any, ...]:
    """Parse atom data without coordinate-precision reporting."""
    syms = block.get('atom_site_type_symbol')
    lbs = block.get('atom_site_label')
    xs = block.get('atom_site_fract_x')
    ys = block.get('atom_site_fract_y')
    zs = block.get('atom_site_fract_z')

    if not isinstance(syms, list) and isinstance(lbs, list):
        # The CIF core dictionary makes _atom_site_type_symbol optional; when it is absent the
        # element is derivable from the label, so infer it rather than refusing the file.
        inferred = [_symbol_from_label(lab) for lab in lbs]
        unresolved = [lab for lab, sym in zip(lbs, inferred) if sym is None]
        if unresolved:
            raise ValueError(
                "CIF block has no _atom_site_type_symbol column and the element could not be "
                f"inferred from _atom_site_label for: {', '.join(unresolved)}"
            )
        warnings.warn(
            "CIF block has no _atom_site_type_symbol column; element symbols inferred from _atom_site_label",
            RuntimeWarning,
            stacklevel=2,
        )
        syms = cast(list[str], inferred)

    missing = [
        f'_{name}'
        for name, column in (
            ('atom_site_type_symbol', syms),
            ('atom_site_label', lbs),
            ('atom_site_fract_x', xs),
            ('atom_site_fract_y', ys),
            ('atom_site_fract_z', zs),
        )
        if not isinstance(column, list)
    ]
    if missing:
        noun = 'column' if len(missing) == 1 else 'columns'
        raise ValueError(f"CIF block is missing required atom-site {noun}: {', '.join(missing)}")
    syms = cast(list[str], syms)
    lbs = cast(list[str], lbs)
    xs = cast(list[str], xs)
    ys = cast(list[str], ys)
    zs = cast(list[str], zs)
    counts = {
        'atom_site_type_symbol': len(syms),
        'atom_site_label': len(lbs),
        'atom_site_fract_x': len(xs),
        'atom_site_fract_y': len(ys),
        'atom_site_fract_z': len(zs),
    }
    if len(set(counts.values())) != 1:
        raise ValueError(
            "CIF atom-site columns have mismatched lengths: " + ", ".join(f"{k}={v}" for k, v in counts.items())
        )

    # Optional occupancy column
    occ_col = block.get('atom_site_occupancy')
    if occ_col is not None:
        occs = []
        occs_exact = []
        occupancy_precisions = []
        for index, t in enumerate(occ_col):
            v, meta = parse_cif_float(t, meta=True)
            occs.append(v)
            occs_exact.append(_prefer_exact_token(block, t, 'httk_atom_site_occupancy_exact', index=index))
            occupancy_precisions.append(combined_precision((meta['precision'], meta['esd'])))
    else:
        occs = None
        occs_exact = None
        occupancy_precisions = None

    symbols = [s.strip() for s in syms]
    labels = [lab.strip() for lab in lbs]

    exact_positions = [
        (
            _prefer_exact_token(block, xi, 'httk_atom_site_fract_x_exact', index=index),
            _prefer_exact_token(block, yi, 'httk_atom_site_fract_y_exact', index=index),
            _prefer_exact_token(block, zi, 'httk_atom_site_fract_z_exact', index=index),
        )
        for index, (xi, yi, zi) in enumerate(zip(xs, ys, zs))
    ]
    positions = [(parse_cif_float(xi), parse_cif_float(yi), parse_cif_float(zi)) for xi, yi, zi in zip(xs, ys, zs)]
    return symbols, labels, positions, exact_positions, occs, occs_exact, occupancy_precisions


def _prefer_exact_token(
    block: Mapping[str, Any], standard: str, companion: str, *, index: int | None = None
) -> str | None:
    """Use an httk exact companion when present, otherwise the standard token."""
    companion_value = block.get(companion)
    if index is not None and isinstance(companion_value, list):
        companion_value = companion_value[index] if index < len(companion_value) else None
    if companion_value is not None:
        token = cif_exact_token(companion_value)
        if token is not None:
            parse_cif_fraction(token)
            return token
    return cif_exact_token(standard)


def _parse_atoms_with_precision(block: Mapping[str, Any]) -> tuple[Any, ...]:
    """Parse atom data and report the coarsest coordinate precision."""
    parsed = _parse_atoms(block)
    xs = block['atom_site_fract_x']
    ys = block['atom_site_fract_y']
    zs = block['atom_site_fract_z']
    companions = tuple(block.get(f'httk_atom_site_fract_{axis}_exact') for axis in 'xyz')
    has_companion = any(value is not None for value in companions)
    claims: list[object] = []
    for index, values in enumerate(zip(xs, ys, zs)):
        for axis, value in enumerate(values):
            companion = companions[axis]
            companion_value = companion[index] if isinstance(companion, list) and index < len(companion) else None
            if companion_value is not None and cif_exact_token(companion_value) is not None:
                continue
            if has_companion and cif_exact_token(value) in {'0', '1'}:
                continue
            entry = parse_cif_float(value, meta=True)[1]
            claims.append(entry['precision'])
            claims.append(entry['esd'])
    return (*parsed, combined_precision(claims))


def _parse_uc_with_precision(block: Mapping[str, Any]) -> tuple[tuple[float | None, ...], Fraction | None]:
    """The six cell parameters, and how precisely the cell was stated.

    The precision is an absolute length, or ``None``, and is taken from the three
    **lengths** only. The angles are deliberately left out: their precision would be in
    degrees rather than the length units everything else here is in, and a right angle in a
    CIF is almost always exact by symmetry rather than a measurement whose last digit means
    anything — folding a "1 degree" precision from a bare ``90`` into a length would be
    nonsense.

    A separate function rather than a flag on :func:`_parse_uc`, because a boolean that
    changes the shape of the return value cannot be narrowed by a type checker.
    """
    tags = _CELL_TAGS
    missing = [tag for tag in tags if block.get(tag) is None]
    if missing:
        raise ValueError(f"CIF block has no unit cell: missing {', '.join('_' + tag for tag in missing)}")

    values = []
    claims: list[object] = []
    for index, tag in enumerate(tags):
        value, entry = parse_cif_float(block[tag], meta=True)
        values.append(value)
        if index < 3:  # lengths only
            token = cif_exact_token(block[tag])
            try:
                exact_integer = token is not None and Fraction(token).denominator == 1
            except (ValueError, ZeroDivisionError):
                exact_integer = False
            if not exact_integer:
                claims.append(entry['precision'])
                claims.append(entry['esd'])
    return tuple(values), combined_precision(claims)


def _basis_from_lengths_angles(
    a: float, b: float, c: float, alpha: float, beta: float, gamma: float
) -> list[list[float]]:
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


class AsuCell(NamedTuple):
    """Hold the parsed cell, atom sites, exact tokens, and precision claims.

    :param basis: Conventional Cartesian basis vectors for the unit cell.
    :param positions: Fractional atom positions.
    :param positions_exact: Original central coordinate tokens for each atom.
    :param occupancies: Parsed atom occupancies, when supplied.
    :param occupancies_exact: Original central occupancy tokens, when supplied.
    :param occupancy_precisions: Coarsest claimed precision for each occupancy, when supplied.
    :param coordinate_precision: Coarsest claimed precision across the coordinates.
    :param basis_precision: Coarsest claimed precision across the cell lengths.
    :param symbols: Atom-site symbols in input order.
    :param labels: Atom-site labels in input order.
    :param equivalent_atoms: One-based identifiers grouping equal atom-site labels.
    """

    basis: list[list[float]]
    positions: list[tuple[float | None, float | None, float | None]]
    positions_exact: list[tuple[str | None, str | None, str | None]]
    occupancies: list[float | None] | None
    occupancies_exact: list[str | None] | None
    occupancy_precisions: list[Fraction | None] | None
    coordinate_precision: Fraction | None
    basis_precision: Fraction | None
    symbols: list[str]
    labels: list[str]
    equivalent_atoms: list[int]


def parse_asu_cell(cifblock: Mapping[str, Any]) -> AsuCell:
    """Parse a CIF block into its asymmetric-unit cell data.

    :param cifblock: Normalized CIF data for one data block.
    :return: The parsed unit cell and atom-site data.
    :raises ValueError: If required cell or atom-site data is missing or invalid.
    """
    parameters, basis_precision = _parse_uc_with_precision(cifblock)
    a, b, c, alpha, beta, gamma = parameters
    if any(value is None for value in parameters):
        raise ValueError("CIF unit-cell parameters cannot be unknown")
    a, b, c, alpha, beta, gamma = cast(tuple[float, float, float, float, float, float], parameters)
    basis = _basis_from_lengths_angles(a, b, c, alpha, beta, gamma)
    (
        symbols,
        labels,
        positions,
        exact_positions,
        occs,
        occs_exact,
        occupancy_precisions,
        coordinate_precision,
    ) = _parse_atoms_with_precision(cifblock)

    # figure out equivalent atoms based on labels
    labels_map = {}
    equivalent_atoms = []
    next_id = 1
    for lab in labels:
        if lab not in labels_map:
            labels_map[lab] = next_id
            next_id += 1
        equivalent_atoms.append(labels_map[lab])

    return AsuCell(
        basis,
        positions,
        exact_positions,
        occs,
        occs_exact,
        occupancy_precisions,
        coordinate_precision,
        basis_precision,
        symbols,
        labels,
        equivalent_atoms,
    )


def parse_structural_modulation(
    cifblock: Mapping[str, Any],
) -> tuple[list[list[Fraction]] | None, int, bool, list[str]]:
    """
    Extract structural superspace modulation information from a standard CIF.

    Returns a tuple ``(structural_q, mod_dim, has_struct_mod, struct_mod_atoms)`` where
        ``structural_q`` is a list of q-vectors or ``None``, ``mod_dim`` is the modulation
        dimension (0 if absent), ``has_struct_mod`` is a bool, and ``struct_mod_atoms`` is a
        sorted list of atom-site labels.

    :param cifblock: Normalized CIF data for one data block.
    :return: Structural q-vectors, modulation dimension, presence flag, and affected labels.
    :raises ValueError: If a structural wave vector contains an unknown component.
    """
    # modulation dimension (0 if absent)
    mod_dim = int(cifblock.get('cell_modulation_dimension', 0))

    # structural_q from cell_wave_vector (only if mod_dim > 0)
    structural_q = None
    qx, qy, qz = (cifblock.get(tag) for tag in CIF_TAGS['structural_q'])
    if qx and qy and qz:
        structural_q = []
        for x, y, z in zip(qx, qy, qz):
            q = [parse_cif_fraction(value) for value in (x, y, z)]
            if any(value is None for value in q):
                raise ValueError("CIF structural wave vector has an unknown component")
            structural_q.append(cast(list[Fraction], q))

    # detect structural Fourier modulations
    has_struct_mod = False
    struct_mod_atoms = set()

    labels = cifblock.get(CIF_TAGS['structural_displacement_label'])
    if labels:
        has_struct_mod = True
        struct_mod_atoms.update(labels)

    labels = cifblock.get(CIF_TAGS['structural_occupancy_label'])
    if labels:
        has_struct_mod = True
        struct_mod_atoms.update(labels)

    return structural_q, mod_dim, has_struct_mod, sorted(struct_mod_atoms)


#: The six unit-cell tags, in the order everything here reports them.
_CELL_TAGS = (
    'cell_length_a',
    'cell_length_b',
    'cell_length_c',
    'cell_angle_alpha',
    'cell_angle_beta',
    'cell_angle_gamma',
)


def _cell_parameter_tokens(cifblock: Mapping[str, Any]) -> tuple[str | None, ...]:
    """The six cell parameters as the text the file wrote, uncertainties stripped.

    The same fidelity argument as ``positions_exact``: ``5.6402`` should be able to become
    the rational 56402/10000, not the binary value of ``float("5.6402")``.
    """
    return tuple(
        _prefer_exact_token(cifblock, tag_value, f'httk_{tag}_exact')
        for tag, tag_value in ((tag, cifblock[tag]) for tag in _CELL_TAGS)
    )


def _first_tag(cifblock: Mapping[str, Any], *names: str) -> Any | None:
    """The first normalized CIF data name present, or ``None``."""
    for name in names:
        value = cifblock.get(name)
        if value is not None:
            return value
    return None


def cifblock_to_asu(cifblock: Mapping[str, Any]) -> dict[str, Any]:
    """Convert one normalized CIF block to a neutral asymmetric-unit payload.

    The payload keeps ``cell_parameters_exact`` and ``positions_exact`` as the central
    numeric text written by the file, and keeps ``symops_xyz`` as the raw operation strings.
    The dual numeric channel gives ``_httk_*_exact`` companion tags precedence over standard
    numeric tags on read, preserving central text after a lossy display write.

    :param cifblock: Normalized CIF data for one data block.
    :return: A neutral mapping containing cell, atom, symmetry, and metadata channels.
    :raises ValueError: If the block lacks required cell, atom-site, or symmetry data.
    """
    asu = parse_asu_cell(cifblock)

    # standard space group symmetry
    symops_xyz = cifblock.get('space_group_symop.operation_xyz')
    if symops_xyz is None:
        # Some readers normalize loop tags without dots, so accept that too.
        symops_xyz = cifblock.get('space_group_symop_operation_xyz')
    if symops_xyz is None:
        # some CIFs use older spelling
        symops_xyz = cifblock.get('symmetry_equiv_pos_as_xyz')

    if symops_xyz is None:
        raise ValueError("CIF block has no symmetry operations")

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
    space_group_nbr = _first_tag(
        cifblock,
        'space_group_it_number',
        'symmetry_space_group_it_number',
        'symmetry_int_tables_number',
    )
    icsd = _first_tag(cifblock, 'database_code_icsd')
    doi = _first_tag(cifblock, 'citation_doi')

    return {
        'format': 'cif',
        'cell_parameters_exact': _cell_parameter_tokens(cifblock),
        'positions_exact': asu.positions_exact,
        'occupancies': asu.occupancies,
        'occupancies_exact': asu.occupancies_exact,
        'occupancy_precisions': asu.occupancy_precisions,
        'symbols': asu.symbols,
        'symops_xyz': tuple(symops_xyz),
        'incomm': incomm,
        'space_group_nbr': space_group_nbr,
        'space_group_name_hm': space_group_name_hm,
        'space_group_name_hall': space_group_name_hall,
        'icsd': icsd,
        'doi': doi,
        'coordinate_precision': asu.coordinate_precision,
        'basis_precision': asu.basis_precision,
        'equivalent_atoms': asu.equivalent_atoms,
        'labels': asu.labels,
    }


def read_cif_asus(source: str | os.PathLike[str] | Iterable[str], *, autocorrect: bool = False) -> dict[str, Any]:
    """Read a CIF and return its asymmetric units as a neutral, tagged payload.

    This is what ``httk.core.load`` returns for a ``.cif`` file: a mapping with
    ``format`` set to ``"cif"``, ``blocks`` holding one asymmetric-unit mapping per data
    block that describes a structure, and ``header`` the file's leading comment. The tag
    lets a consumer dispatch on the file type without knowing which reader
    produced the payload. When ``autocorrect=True``, the top-level payload also contains
    ``autocorrect=True`` so downstream adapters can apply compatible repairs.

    Loading never fails on account of a block that is not a structure. CIF is a
    general-purpose format and a file may hold bibliographic entries, powder patterns, or
    an incomplete draft alongside — or instead of — anything crystallographic. Blocks
    without atom sites are simply not structures and are passed over; blocks that have
    atom sites but cannot be interpreted are collected in ``unparsed``, each with the
    reason, so that nothing is dropped silently and the failure surfaces when a structure
    is actually asked for.

    The mapping stays neutral — plain lists, strings, exact numeric text channels, and raw
    symmetry-operation strings, with no domain objects — so *httk-io* need not know about
    *httk-atomistic*. Turning it into a structure is
    ``httk.core.load`` (which returns an ``ASUStructure`` when atomistic support
    is installed).

    :param source: A filename, open text stream, or iterable of CIF lines.
    :param autocorrect: Drop malformed auxiliary loops and warn about each repair.
    :return: A neutral CIF payload containing structural blocks, unparsed block reasons, and the header.
    :raises ValueError: If the CIF stream contains malformed data that prevents parsing.
    """
    cifblocks, header = read_cif(source, allow_cif2=False, autocorrect=autocorrect)

    blocks = []
    unparsed = []
    for name, cifblock in cifblocks:
        if 'atom_site_label' not in cifblock:
            continue
        try:
            blocks.append(cifblock_to_asu(cifblock))
        except Exception as error:
            unparsed.append({'block': name, 'reason': f'{type(error).__name__}: {error}'})

    payload: dict[str, Any] = {'format': 'cif', 'blocks': blocks, 'unparsed': unparsed, 'header': header}
    if autocorrect:
        payload['autocorrect'] = True
    return payload


def single_asu_from_cif_file(source: str | os.PathLike[str] | Iterable[str]) -> dict[str, Any]:
    """Return the first structural CIF block from :func:`read_cif_asus`.

    :param source: A filename, open text stream, or iterable of CIF lines.
    :return: The first parsed asymmetric-unit mapping.
    :raises ValueError: If no structural block is available.
    """
    payload = read_cif_asus(source)
    if payload['blocks']:
        return payload['blocks'][0]
    if payload['unparsed']:
        raise ValueError(payload['unparsed'][0]['reason'])
    raise ValueError("No structural block found in CIF")
