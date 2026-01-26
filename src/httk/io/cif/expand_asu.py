import itertools
import numpy as np
from .cif_parser import single_asu_from_cif_file
from .mcif_parser import single_mag_asu_from_mcif_file

def _wrap01(v):
    v = np.asarray(v, float)
    return v - np.floor(v)

def _snap_to_grid(v, grid_den, tol=1e-8):
    """
    Snap fractional coords to nearest multiple of 1/grid_den, wrap to [0,1).
    """
    v = _wrap01(v)
    idx = np.floor(v * grid_den + 0.5).astype(np.int64)
    idx %= grid_den
    return idx / grid_den

def pos_key(v, grid_den, tol=1e-8):
    """Hashable key for a position after snapping."""
    v_snap = _snap_to_grid(np.asarray(v, float), grid_den, tol)
    return tuple(int(round(x * grid_den)) for x in v_snap)

def _wrap_neighbor_key(k, delta, grid_dens):
    """Return neighbor key with periodic wrap-around."""
    return tuple((ki + di) % grid_dens for ki, di in zip(k, delta))

def find_adjacent_key_wrap(k, pos_map, grid_dens):
    """
    Look for k in pos_map; if missing, check the 26 neighbors
    (Chebyshev distance 1) with periodic boundaries.
    Return (found_key or None).
    """
    if k in pos_map:
        return k

    for delta in itertools.product((-1, 0, 1), repeat=3):
        if delta == (0, 0, 0):
            continue
        k2 = _wrap_neighbor_key(k, delta, grid_dens)
        if k2 in pos_map:
            return k2

    return None

def apply_op_frac(R, t, f):
    """
    Apply a (fractional) symmetry operation and wrap to [0,1).
    Works with int-like R and float-like t/f.
    """
    R = np.asarray(R, int)
    t = np.asarray(t, float)
    f = np.asarray(f, float)
    return _wrap01(R @ f + t)

# ---------- Species helpers ----------

def species_to_numbers(species):
    unique_species = sorted(set(species))
    species_to_number = {sp: i + 1 for i, sp in enumerate(unique_species)}
    numbers_to_species = {i + 1: sp for i, sp in enumerate(unique_species)}
    return [species_to_number[sp] for sp in species], numbers_to_species

def asu_data_to_numbers_by_labels(asu_data):
    """
    Optional helper mirroring your magnetic 'separate_noneq_atoms' behavior:
    give each distinct *label* a unique number (even if same element symbol).
    """
    labels_map = {}
    numbers_to_species = {}
    numbers = []
    next_id = 1
    for sym, lab in zip(asu_data["symbols"], asu_data["labels"]):
        if lab not in labels_map:
            labels_map[lab] = next_id
            numbers_to_species[next_id] = sym
            next_id += 1
        numbers.append(labels_map[lab])
    return numbers, numbers_to_species

# ---------- Generalized expander ----------

def expand_asu(
    positions_frac,
    species,
    ops,
    *,
    grid_dens=16384,
    tol=1e-6,
    moments=None,
    apply_op_moment_fn=None,
    spin_basis=None,
    basis=None,
    mag_grid_dens=16384,
):
    """
    Generalized ASU expansion.

    Nonmag usage:
        positions_full, species_full = expand_asu_general(pos, species, ops)

    Mag usage:
        positions_full, species_full, moments_full, mom_classes, symops_classes = expand_asu_general(
            pos, species, ops, moments=moments, apply_op_moment_fn=apply_op_moment,
            spin_basis=spin_basis, basis=basis, mag_grid_dens=mag_grid_dens
        )

    Notes
    -----
    - For nonmag, `ops` are (R, t).
    - For mag,  `ops` are (R, t, time_flag).
    - Dedup logic matches your original (snap keys + neighbor wrap).
    """

    positions_frac = list(positions_frac)
    species = list(species)

    if moments is not None:
        moments = list(moments)
        if len(positions_frac) != len(species) or len(species) != len(moments):
            raise ValueError("len(positions_frac), len(species), len(moments) must match")
        if apply_op_moment_fn is None:
            raise ValueError("moments provided but apply_op_moment_fn is None")
        if spin_basis is None:
            raise ValueError("moments provided but spin_basis is None")
        if basis is None:
            raise ValueError("moments provided but basis is None")

    # For nonmag
    seen_by_pos = set()
    pos_full = []
    species_full = []

    # For mag (kept compatible with your original outputs)
    mom_full = []
    mom_classes = []
    mom_dirs = None
    ops_map = None

    if moments is not None:
        # replicate your "moment direction classes" machinery
        if spin_basis in ["cartesian", "crystal"]:
            zero_spin = (0, 0, 0)
        elif spin_basis in ["collinear"]:
            zero_spin = 0
        else:
            raise Exception("Unexpected spin_basis: " + str(spin_basis))

        mom_dirs = [zero_spin]
        ops_map = [{} for _ in range(len(ops))]

        def get_snapped_magdir(v, mag_grid_dens, spin_basis, *, tol=1e-8, int_mult_of_grid=False):
            # directly lifted from your snippet (kept logic)
            if spin_basis == "collinear":
                if abs(v) < tol:
                    return 0, 0
                return 1, 1 if v > 0 else -1

            v = np.asarray(v, dtype=float)
            n = np.linalg.norm(v)
            if n < tol:
                return (0, 0, 0), 0
            magdir = v / n
            k = np.trunc(magdir * mag_grid_dens).astype(int)

            if not np.any(k):
                return (0, 0, 0), 0

            first_idx = np.flatnonzero(k)[0]
            sig = 1 if k[first_idx] > 0 else -1
            if sig < 0:
                k = -k

            if int_mult_of_grid:
                return tuple(k.tolist()), sig

            return k / mag_grid_dens, sig

        # map from snapped-position key -> (moment, moment_class)
        seen_by_pos_mag = {}

    for seed_idx, (seed_pos, specie) in enumerate(zip(positions_frac, species)):
        f = np.asarray(seed_pos, float)
        if moments is not None:
            m0 = moments[seed_idx]
            m0dir, m0dir_sign = get_snapped_magdir(m0, mag_grid_dens, spin_basis, int_mult_of_grid=True)
            if m0dir not in mom_dirs:
                mom_dirs.append(m0dir)
            m0class = mom_dirs.index(m0dir) * m0dir_sign

        for ops_idx, op in enumerate(ops):
            if moments is None:
                # nonmag: op=(R,t)
                R, t = op
                g = apply_op_frac(R, t, f)
            else:
                # mag: op=(R,t,time_flag)
                R, t, time_flag = op
                g = apply_op_frac(R, t, f)
                m = apply_op_moment_fn(R, time_flag, m0, basis, spin_basis)

                k = pos_key(g, grid_dens, tol=1e-8)
                kf = find_adjacent_key_wrap(k, seen_by_pos_mag, grid_dens)

                mdir, mdir_sign = get_snapped_magdir(m, mag_grid_dens, spin_basis, int_mult_of_grid=True)
                if mdir not in mom_dirs:
                    mom_dirs.append(mdir)
                mclass = mom_dirs.index(mdir) * mdir_sign

                # op spin mapping consistency (kept)
                if mdir != (0, 0, 0) and m0class != mclass:
                    if m0class in ops_map[ops_idx]:
                        if ops_map[ops_idx][m0class] != mclass:
                            raise Exception(
                                "Inconsistent op spin mapping: "
                                + str(m0class) + "->" + str(ops_map[ops_idx][m0class])
                                + " vs. " + str(mclass) + " for: " + str(op)
                            )
                    else:
                        ops_map[ops_idx][m0class] = mclass

                if kf is not None and kf in seen_by_pos_mag:
                    m_prev, m_prev_class = seen_by_pos_mag[kf]
                    if not np.allclose(m, m_prev, atol=1e-6, rtol=0):
                        if spin_basis != "collinear":
                            raise Exception(
                                "Internally inconsistent positions and symmetry ops: "
                                + str(op) + " transforms " + str(f) + ":" + str(m0)
                                + " into " + str(g) + ":" + str(m)
                                + " but site is: " + str(m_prev)
                            )
                        elif abs(m) - abs(m_prev) < 1e-6:
                            m = m_prev
                            mclass = m_prev_class
                        else:
                            raise Exception(
                                "Internally inconsistent positions and symmetry ops: "
                                + str(op) + " transforms " + str(f) + ":" + str(m0)
                                + " into " + str(g) + ":" + str(m)
                                + " but site is: " + str(m_prev)
                            )
                    continue

                seen_by_pos_mag[k] = (m, mclass)
                pos_full.append(g)
                species_full.append(specie)
                mom_full.append(tuple(m) if spin_basis != "collinear" else m)
                mom_classes.append(mclass)
                continue

            # ---- Nonmag branch dedup ----
            k = pos_key(g, grid_dens, tol=1e-8)
            kf = find_adjacent_key_wrap(k, seen_by_pos, grid_dens)
            if kf is not None and kf in seen_by_pos:
                continue

            seen_by_pos.add(k)
            pos_full.append(g)
            species_full.append(specie)

    # Sort for reproducibility (same style as you used)
    order = sorted(
        range(len(pos_full)),
        key=lambda i: (round(pos_full[i][2], 6), round(pos_full[i][1], 6), round(pos_full[i][0], 6)),
    )
    pos_full = [pos_full[i] for i in order]
    species_full = [species_full[i] for i in order]

    if moments is None:
        return pos_full, species_full

    mom_full = [mom_full[i] for i in order]
    mom_classes = [mom_classes[i] for i in order]
    symops_classes = tuple(tuple(sorted(m.items())) for m in ops_map)
    return pos_full, species_full, mom_full, mom_classes, symops_classes


# ---------- Nonmag CIF convenience ----------

def _apply_op_moment_crystal(R, time_flag, m):
     """
     Magnetic moment is an axial vector: m' = time * det(R) * R * m
     R is the same 3x3 integer matrix used on fractional coords (active rotation in the crystal basis).
     """
     R = np.asarray(R, int)
     m = np.asarray(m, float)
     detR = _det_int3(R)  # ±1
     if time_flag != 0:
         return float(time_flag) * detR * (R @ m)
     else:
         return detR * (R @ m)

def _apply_op_moment_cartesian(R, time_reversal, m, lattice):
    """
    R: 3x3 int rotation in fractional basis (spglib's rotations[k])
    time_reversal: 0 or 1 (spglib's time_reversals[k])
    m: 3-vector spin (Cartesian)
    lattice: 3x3 lattice with row vectors [a; b; c] in Cartesian (spglib style)
    """
    R = np.asarray(R, int)
    m = np.asarray(m, float)

    # Build A with lattice vectors as columns
    A = np.asarray(lattice, float).T

    # Cartesian rotation equivalent to W
    R_cart = A @ R @ np.linalg.inv(A)

    # Axial-vector extra sign for improper ops
    detR = int(round(np.linalg.det(R)))  # ±1 for spglib rotations
    s = (detR * R_cart) @ m

    # Time reversal flips spin
    #if time_reversal != 0:
    #    s = float(time_reversal)*s
    # Spglib convention
    if time_reversal == 1:
        s = -1.0*s
    elif time_reversal != 0:
        raise Exception("Inconsitency error: unexpected time reversal flag:"+str(time_reversal))

    # Optional: clean tiny numerical noise
    s[np.isclose(s, 0.0, atol=1e-12)] = 0.0
    return s

def _apply_op_moment_collinear(R, time_reversal, m0):
    """
    Simplified collinear spin (basis-0 scalar) transform.
    Only tracks 'up'/'down' sign; ignores axis/lattice.

    Parameters
    ----------
    R : (3,3) int-like
        Rotation in fractional basis (spglib rotations[k]).
    time_reversal : int
        0 or 1 (spglib time_reversals[k]).
    m0 : float
        Scalar moment (+ for up, - for down).

    Returns
    -------
    m0_prime : float
        Transformed scalar moment (sign-preserved or flipped).

    Rule
    ----
    q = (-1)^tau * det(R).
    If q = +1 -> keep sign; if q = -1 → flip sign.

    Notes
    -----
    This intentionally ignores whether R rotates the collinear axis.
    For precise axis-aware behavior, use the basis-1 version.
    """
    if time_reversal not in (0, 1):
        raise ValueError(f"Unexpected time_reversal flag: {time_reversal!r}")

    R = np.asarray(R, int)
    detR = int(round(np.linalg.det(R)))
    #if detR not in (+1, -1):
    #    raise ValueError(f"det(R) must be +/-1, got {detR}")

    flips = (time_reversal + (1 if detR == -1 else 0)) % 2
    return -m0 if flips else m0

def apply_op_moment(R, time_flag, m, lattice, spin_basis):
    if spin_basis == "collinear":
        return _apply_op_moment_collinear(R, time_flag, m)
    if spin_basis == "cartesian":
        return _apply_op_moment_cartesian(R, time_flag, m, lattice)
    if spin_basis == "crystal":
        return _apply_op_moment_crystal(R, time_flag, m)
    raise Exception("Unrecognized spin_basis: "+str(spin_basis))

def mag_asu_data_to_numbers(mag_asu_data):
    """
    Numbering that separates non-equivalent atoms (labels), mapping each label to a unique id,
    while still storing element symbol as numbers_to_species[id].
    This is the corrected intent of your original helper.
    """
    label_to_id = {}
    numbers_to_species = {}
    numbers = []
    next_id = 1
    for sym, lab in zip(mag_asu_data["symbols"], mag_asu_data["labels"]):
        if lab not in label_to_id:
            label_to_id[lab] = next_id
            numbers_to_species[next_id] = sym
            next_id += 1
        numbers.append(label_to_id[lab])
    return numbers, numbers_to_species


# -------------

def cif_to_struct(fs, *, separate_noneq_atoms=False, tol=1e-6):
    """
    Read a CIF file, expand the asymmetric unit, and return the result

    Parameters
    ----------
    fs : filename or file-like (as accepted by your read_cif layer)
    separate_noneq_atoms : bool
        If True, assign unique numbers by CIF labels (even if same element).
        If False, assign by element symbol.
    tol : float
        Used in expansion consistency (position hashing uses fixed 1e-8 like your original).

    Returns
    -------
    dict with data
    """

    asu_data = single_asu_from_cif_file(fs)

    # Use the resolution returned by _parse_atoms(...) from cif_parser layer
    grid_dens = int(1.0 / asu_data["resolution"]) if asu_data.get("resolution") not in (None, 0.0) else 16384

    # Expand: note CIF symops are (R,t) without time flag
    positions_full, species_full = expand_asu(
        asu_data["positions"],
        asu_data["symbols"] if not separate_noneq_atoms else asu_data["labels"],
        asu_data["symops"],
        grid_dens=grid_dens,
        tol=tol,
    )

    if separate_noneq_atoms:
        # each label is its own "species number", but map back to element symbol
        # numbers_to_species: number -> element symbol
        numbers, numbers_to_species = asu_data_to_numbers_by_labels(asu_data)
        # expand produced labels_full; convert those to numbers via label->id
        label_to_id = {lab: i for i, lab in enumerate(sorted(set(asu_data["labels"])), start=1)}
        numbers_full = [label_to_id[l] for l in species_full]
        species_meta = [numbers_to_species[i] for i in numbers_full]
    else:
        numbers_full, numbers_to_species = species_to_numbers(species_full)
        species_meta = species_full

    meta = {
        "species": species_meta,
        "numbers_to_species": numbers_to_species,
        "grid_dens": grid_dens,
    }
    asu_data.update(meta)

    asu_data["positions_full"] = positions_full
    asu_data["numbers_full"] = numbers_full
    asu_data["species_full"] = species_full

    return asu_data


def expand_mag_asu(fs, *, separate_noneq_atoms=False, error_on_nonmag=True, tol=1e-6, mag_grid_dens=16384):
    """
    Read an mCIF file, expand the asymmetric unit, and return the result

    Parameters
    ----------
    fs : filename or file-like
        As accepted by your read_cif() stack.
    separate_noneq_atoms : bool
        If True, assign unique numbers by CIF labels (even if same element symbol).
        If False, assign by element symbol.
    error_on_nonmag : bool
        If True, raise if no magnetic moments can be extracted.
    tol : float
        Expansion consistency tolerance (matches your original usage).
    mag_grid_dens : int
        Density used when snapping magnetization directions/classes (passed through the generalized expander).

    Returns
    -------
    dict with data
    """
    mag_asu_data = single_mag_asu_from_mcif_file(fs, error_on_nonmag=error_on_nonmag)

    # Same idea as your original:
    grid_dens = int(1.0 / mag_asu_data["resolution"]) if mag_asu_data.get("resolution") not in (None, 0.0) else 16384

    if separate_noneq_atoms:
        # Use label-numbering for the ASU seeds, so equivalent atoms can be tracked separately if desired
        asu_numbers, numbers_to_species = mag_asu_data_to_numbers(mag_asu_data)

        positions_full, numbers_full, moments_full, momclasses, symops_classes = expand_asu(
            mag_asu_data["positions"],
            asu_numbers,
            mag_asu_data["symops"],
            grid_dens=grid_dens,
            mag_grid_dens=mag_grid_dens,
            tol=tol,
            moments=mag_asu_data["magmoms"],
            apply_op_moment_fn=apply_op_moment,
            spin_basis=mag_asu_data["spin_basis"],
            basis=mag_asu_data["basis"],
        )

        species_full = [numbers_to_species[i] for i in numbers_full]
    else:
        positions_full, species_full, moments_full, momclasses, symops_classes = expand_asu(
            mag_asu_data["positions"],
            mag_asu_data["symbols"],
            mag_asu_data["symops"],
            grid_dens=grid_dens,
            mag_grid_dens=mag_grid_dens,
            tol=tol,
            moments=mag_asu_data["magmoms"],
            apply_op_moment_fn=apply_op_moment,
            spin_basis=mag_asu_data["spin_basis"],
            basis=mag_asu_data["basis"],
        )

        numbers_full, numbers_to_species = species_to_numbers(species_full)

    # Add metadata (same spirit as your old function)
    meta = {
        "species": species_full,
        "numbers_to_species": numbers_to_species,
        "grid_dens": grid_dens,
        "momclasses": momclasses,
        "symops_classes": symops_classes,
        "mag_grid_dens": mag_grid_dens,
    }
    mag_asu_data.update(meta)

    mag_asu_data["positions_full"] = positions_full
    mag_asu_data["numbers_full"] = numbers_full
    mag_asu_data["moments_full"] = moments_full
    mag_asu_data["species_full"] = species_full

    return mag_asu_data
