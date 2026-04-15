"""
BBGKY hierarchy propagator for the quantum spin-boson / vibronic coupling model.

Truncation schemes
------------------
L1 (default)
    Include all (m1, m2) with  sum(m1) + sum(m2) < tier.
    N grows as  O(tier^{2·ndof})  — use for small ndof or small tier.

Per-mode  (mode_dims keyword)
    Include all (m1, m2) with  m1[k] + m2[k] <= mode_dims[k]  for every k.
    N = prod_k  (mode_dims[k]+1)(mode_dims[k]+2)/2
    Better than L1 when coupling strengths differ strongly across modes
    (assign smaller mode_dims to weakly-coupled or high-frequency modes).

Filtering + active-set propagation
    Call filter_ados() every M steps to zero out negligibly small ADOs,
    then pass the returned active mask to prop().  Only ADOs that are
    non-zero OR adjacent to a non-zero ADO are propagated, giving a
    speedup proportional to the sparsity of the hierarchy at runtime.

MPS / Tensor-Train note (not implemented here)
    For ndof >= 5, represent the full ADO tensor as a Matrix Product State
    (bond dimension chi).  Cost scales as O(ndof * chi^3) instead of
    O(tier^{2*ndof}).  See: Borrelli & Gelin, J. Chem. Phys. 150, 234102
    (2019); Ren et al., J. Phys. Chem. Lett. 12, 2262 (2021).
"""

import numpy as np
from itertools import combinations, product as iproduct


# ──────────────────────────────────────────────────────────────────────────────
# ADO hierarchy construction
# ──────────────────────────────────────────────────────────────────────────────

def _multiindex_list(ndof: int, total: int):
    """All ndof-tuples of non-negative integers that sum to `total`."""
    if ndof == 1:
        return [(total,)]
    result = []
    for combo in combinations(range(ndof + total - 1), ndof - 1):
        vec = [combo[0]]
        for j in range(1, ndof - 1):
            vec.append(combo[j] - combo[j - 1] - 1)
        vec.append(ndof + total - 2 - combo[-1])
        result.append(tuple(vec))
    return result


def build_ado_structure(ndof: int, tier: int = None, mode_dims=None):
    """
    Enumerate ADO multi-indices and precompute neighbour lookup tables.

    Exactly one of `tier` or `mode_dims` must be supplied.

    Parameters
    ----------
    ndof      : number of bath modes
    tier      : L1 truncation level — include ADOs with sum(m1)+sum(m2) < tier
    mode_dims : per-mode max combined occupation, array-like of length ndof.
                ADO (m1, m2) is included iff m1[k]+m2[k] <= mode_dims[k] for all k.
                Example for ndof=3: mode_dims=[5, 5, 4]

    Returns
    -------
    N          : total number of ADOs
    up1, up2   : (N, ndof) int arrays — neighbours at (m1+ek, m2) and (m1, m2+ek)
    dn1, dn2   : (N, ndof) int arrays — neighbours at (m1-ek, m2) and (m1, m2-ek)
    m1_arr     : (N, ndof) float array of m1 vectors
    m2_arr     : (N, ndof) float array of m2 vectors

    Sentinel index N maps to the zero-padded extra slice that must be appended
    to rhot before calling _compute_rhs.
    """
    if (tier is None) == (mode_dims is None):
        raise ValueError("Provide exactly one of `tier` or `mode_dims`.")

    if mode_dims is not None:
        # ── per-mode constraint ────────────────────────────────────────────
        mode_dims = np.asarray(mode_dims, dtype=int)
        # For mode k: all (m1k, m2k) with m1k + m2k <= mode_dims[k]
        per_mode = [
            [(a, b) for a in range(d + 1) for b in range(d + 1 - a)]
            for d in mode_dims
        ]
        ado_list = [
            (tuple(p[0] for p in combo), tuple(p[1] for p in combo))
            for combo in iproduct(*per_mode)
        ]
    else:
        # ── L1 constraint (original) ──────────────────────────────────────
        ado_list = []
        for k in range(tier):
            for l in range(k + 1):
                for m1 in _multiindex_list(ndof, k - l):
                    for m2 in _multiindex_list(ndof, l):
                        ado_list.append((m1, m2))

    N    = len(ado_list)
    ZERO = N  # sentinel: out-of-hierarchy → zero-padded slice

    index_map = {key: i for i, key in enumerate(ado_list)}

    up1    = np.full((N, ndof), ZERO, dtype=np.intp)
    up2    = np.full((N, ndof), ZERO, dtype=np.intp)
    dn1    = np.full((N, ndof), ZERO, dtype=np.intp)
    dn2    = np.full((N, ndof), ZERO, dtype=np.intp)
    m1_arr = np.zeros((N, ndof), dtype=float)
    m2_arr = np.zeros((N, ndof), dtype=float)

    for i, (m1, m2) in enumerate(ado_list):
        m1_arr[i] = m1
        m2_arr[i] = m2
        for k in range(ndof):
            nm1p = m1[:k] + (m1[k] + 1,) + m1[k + 1:]
            nm2p = m2[:k] + (m2[k] + 1,) + m2[k + 1:]
            up1[i, k] = index_map.get((nm1p, m2), ZERO)
            up2[i, k] = index_map.get((m1, nm2p), ZERO)
            if m1[k] > 0:
                nm1n = m1[:k] + (m1[k] - 1,) + m1[k + 1:]
                dn1[i, k] = index_map.get((nm1n, m2), ZERO)
            if m2[k] > 0:
                nm2n = m2[:k] + (m2[k] - 1,) + m2[k + 1:]
                dn2[i, k] = index_map.get((m1, nm2n), ZERO)

    return N, up1, up2, dn1, dn2, m1_arr, m2_arr


# ──────────────────────────────────────────────────────────────────────────────
# Importance filtering
# ──────────────────────────────────────────────────────────────────────────────

def filter_ados(rhot, threshold=1e-12):
    """
    Zero out auxiliary density operators whose Frobenius norm is below
    `threshold`.  The zeroth ADO (physical reduced density matrix) is
    always protected.

    Parameters
    ----------
    rhot      : (S, S, N) complex array, modified in-place
    threshold : Frobenius norm cut-off

    Returns
    -------
    active : (N,) boolean array — True for ADOs that were kept
    """
    norms  = np.linalg.norm(rhot.reshape(-1, rhot.shape[2]), axis=0)  # (N,)
    active = norms >= threshold
    active[0] = True                    # protect physical density matrix
    rhot[:, :, ~active] = 0.0
    return active


def _build_work_idx(active, up1, up2, dn1, dn2):
    """
    Return indices of ADOs that need RHS evaluation.

    An ADO needs evaluating if it is active (non-negligible) OR if any of
    its hierarchy neighbours is active (it may receive a non-zero coupling
    contribution in the next step).
    """
    N          = len(active)
    active_ext = np.append(active, False)   # sentinel at position N → False

    work = active.copy()
    for k in range(up1.shape[1]):
        work |= active_ext[up1[:, k]]
        work |= active_ext[up2[:, k]]
        work |= active_ext[dn1[:, k]]
        work |= active_ext[dn2[:, k]]

    return np.where(work)[0]


# ──────────────────────────────────────────────────────────────────────────────
# Propagator
# ──────────────────────────────────────────────────────────────────────────────

def _compute_rhs(rhot_p, Hs, Qs, gamma, up1, up2, dn1, dn2,
                 m1_arr, m2_arr, coeff, ndof, work_idx):
    """
    Compute drho/dt for the ADOs specified by work_idx.

    Parameters
    ----------
    rhot_p   : (S, S, N+1) — last slice is the zero-sentinel
    work_idx : 1-D int array of ADO indices to evaluate

    Returns
    -------
    rhs : (S, S, len(work_idx)) complex array
    """
    rho = rhot_p[:, :, work_idx].transpose(2, 0, 1)   # (W, S, S)

    rhs  = -1j * (Hs @ rho - rho @ Hs)
    rhs += gamma[work_idx, np.newaxis, np.newaxis] * rho

    for k in range(ndof):
        rho_u1 = rhot_p[:, :, up1[work_idx, k]].transpose(2, 0, 1)
        rho_u2 = rhot_p[:, :, up2[work_idx, k]].transpose(2, 0, 1)
        rho_d1 = rhot_p[:, :, dn1[work_idx, k]].transpose(2, 0, 1)
        rho_d2 = rhot_p[:, :, dn2[work_idx, k]].transpose(2, 0, 1)

        rho_up = rho_u1 + rho_u2
        m1k    = m1_arr[work_idx, k, np.newaxis, np.newaxis]
        m2k    = m2_arr[work_idx, k, np.newaxis, np.newaxis]
        Q  = Qs[k]
        c2 = coeff[k] ** 2
        rhs += -1j * (Q @ rho_up - rho_up @ Q)
        rhs += (-0.5j * c2) * (m1k * (Q @ rho_d1) - m2k * (rho_d2 @ Q))

    return rhs.transpose(1, 2, 0)   # (S, S, W)


def _append_zero_sentinel(rhot):
    """Append one zero slice used by out-of-hierarchy neighbour lookups."""
    return np.concatenate(
        [rhot, np.zeros((rhot.shape[0], rhot.shape[1], 1), dtype=complex)],
        axis=2,
    )


def _stage_state(rhot, work_idx, delta):
    """Return a stage state where only `work_idx` slices are incremented."""
    stage = rhot.copy()
    stage[:, :, work_idx] += delta
    return stage


def prop(rhot, parameters, active=None):
    """
    Advance ADOs by one 4th-order Taylor step (RK4 for this linear ODE).

    Parameters
    ----------
    rhot       : (S, S, N) complex array, updated in-place and returned
    parameters : parameter object — must carry Hs, Qs, c, dt, ndof, nmod,
                 up1, up2, dn1, dn2, m1_arr, m2_arr, gamma
    active     : (N,) boolean mask returned by filter_ados().
                 When provided, only active ADOs and their immediate
                 hierarchy neighbours are propagated (sparse update).
                 Pass None (default) to propagate all N ADOs.

    Returns
    -------
    rhot : same array, updated in-place
    """
    Hs    = parameters.Hs
    Qs    = parameters.Qs
    ndof  = parameters.ndof
    coeff = parameters.c
    dt    = parameters.dt
    up1, up2 = parameters.up1, parameters.up2
    dn1, dn2 = parameters.dn1, parameters.dn2
    m1_arr   = parameters.m1_arr
    m2_arr   = parameters.m2_arr
    gamma    = parameters.gamma

    N = rhot.shape[2]

    # Build the work set ───────────────────────────────────────────────────
    if active is None:
        work_idx = np.arange(N, dtype=np.intp)
    else:
        work_idx = _build_work_idx(active, up1, up2, dn1, dn2)

    coeff = np.asarray(coeff, dtype=float)
    if coeff.ndim == 2:
        if coeff.shape == (ndof, ndof) and np.allclose(coeff, np.diag(np.diag(coeff))):
            coeff = np.diag(coeff)
        elif coeff.shape == (1, ndof):
            coeff = coeff[0]
        elif coeff.shape == (ndof, 1):
            coeff = coeff[:, 0]
        else:
            raise ValueError(
                "BBGKY expects one scalar coupling per bath mode; "
                "provide `c` as a length-ndof vector or diagonal matrix."
            )
    if coeff.shape != (ndof,):
        raise ValueError(f"`c` must have shape ({ndof},), got {coeff.shape}.")

    Qs = np.asarray(Qs, dtype=complex)
    if Qs.shape != (ndof, Hs.shape[0], Hs.shape[1]):
        raise ValueError(
            "BBGKY expects one system operator per bath mode; "
            f"`Qs` must have shape ({ndof}, {Hs.shape[0]}, {Hs.shape[1]}), got {Qs.shape}."
        )

    rhot_p = _append_zero_sentinel(rhot)
    k1 = _compute_rhs(rhot_p, Hs, Qs, gamma, up1, up2, dn1, dn2,
                      m1_arr, m2_arr, coeff, ndof, work_idx)

    rhot_2 = _stage_state(rhot, work_idx, 0.5 * dt * k1)
    rhot_2p = _append_zero_sentinel(rhot_2)
    k2 = _compute_rhs(rhot_2p, Hs, Qs, gamma, up1, up2, dn1, dn2,
                      m1_arr, m2_arr, coeff, ndof, work_idx)

    rhot_3 = _stage_state(rhot, work_idx, 0.5 * dt * k2)
    rhot_3p = _append_zero_sentinel(rhot_3)
    k3 = _compute_rhs(rhot_3p, Hs, Qs, gamma, up1, up2, dn1, dn2,
                      m1_arr, m2_arr, coeff, ndof, work_idx)

    rhot_4 = _stage_state(rhot, work_idx, dt * k3)
    rhot_4p = _append_zero_sentinel(rhot_4)
    k4 = _compute_rhs(rhot_4p, Hs, Qs, gamma, up1, up2, dn1, dn2,
                      m1_arr, m2_arr, coeff, ndof, work_idx)

    rhot[:, :, work_idx] += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return rhot
