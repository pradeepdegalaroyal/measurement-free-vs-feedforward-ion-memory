"""Sweeps and statistics built on the exact schemes."""
from __future__ import annotations

from dataclasses import asdict

import numpy as np

from .schemes import Params, measurement_based_fail, measurement_free_fail

TAU_LO, TAU_HI = 0.5, 1.0e4


def tau_star(prm: Params, mf_fail: float | None = None) -> float:
    """Measurement time (gate-times) at which feed-forward and measurement-free tie.

    The feed-forward failure grows monotonically with ``tau_m`` while the measurement-free
    failure does not depend on it, so bisection on ``log tau_m`` finds the crossing.
    Returns ``TAU_LO`` if measurement-free wins even for the fastest measurement and
    ``TAU_HI`` if feed-forward still wins at the largest tested time.
    """
    target = measurement_free_fail(prm) if mf_fail is None else mf_fail
    if measurement_based_fail(prm.with_(tau_m=TAU_LO)) >= target:
        return TAU_LO
    if measurement_based_fail(prm.with_(tau_m=TAU_HI)) <= target:
        return TAU_HI
    lo, hi = np.log(TAU_LO), np.log(TAU_HI)
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if measurement_based_fail(prm.with_(tau_m=float(np.exp(mid)))) < target:
            lo = mid
        else:
            hi = mid
    return float(np.exp(0.5 * (lo + hi)))


def phase_map(base: Params, p2_grid: np.ndarray, tau_grid: np.ndarray) -> dict[str, np.ndarray]:
    """Logical failure of both schemes on a ``(p2, tau_m)`` grid.

    Returns arrays of shape ``(len(tau_grid), len(p2_grid))``.
    """
    mb = np.zeros((len(tau_grid), len(p2_grid)))
    mf = np.zeros_like(mb)
    for j, p2 in enumerate(p2_grid):
        prm = base.with_(p2=float(p2))
        mf[:, j] = measurement_free_fail(prm)
        for i, tau in enumerate(tau_grid):
            mb[i, j] = measurement_based_fail(prm.with_(tau_m=float(tau)))
    return {"mb": mb, "mf": mf}


def shots_to_resolve(p_a: float, p_b: float, z_alpha: float = 1.96, z_beta: float = 0.8416) -> float:
    """Shots per arm for a two-sided two-proportion test (95 % level, 80 % power)."""
    delta = abs(p_a - p_b)
    if delta == 0.0:
        return float("inf")
    pbar = 0.5 * (p_a + p_b)
    num = z_alpha * np.sqrt(2 * pbar * (1 - pbar)) + z_beta * np.sqrt(
        p_a * (1 - p_a) + p_b * (1 - p_b))
    return float((num / delta) ** 2)


def finite_shot_difference(p_mb: float, p_mf: float, shots: int, seeds: list[int]) -> dict[str, np.ndarray]:
    """Binomial emulation of an experiment: estimated ``P_MB - P_MF`` for each seed."""
    est, z = [], []
    for s in seeds:
        rng = np.random.default_rng(1000 + s)
        a = rng.binomial(shots, p_mb) / shots
        b = rng.binomial(shots, p_mf) / shots
        se = np.sqrt(a * (1 - a) / shots + b * (1 - b) / shots)
        est.append(a - b)
        z.append((a - b) / se if se > 0 else 0.0)
    return {"diff": np.array(est), "z": np.array(z)}


def calibration_regret(base: Params, p2_cells: np.ndarray, tau_cells: np.ndarray,
                       sigma: float, seed: int) -> dict[str, np.ndarray]:
    """Scheme-choice regret when the noise model used to decide is mis-characterised.

    Each cell draws log-normal multiplicative errors (spread ``sigma``) on ``p2``,
    ``eps0``, ``p_meas`` and ``tau_m``; the scheme with the smaller *believed* failure
    is chosen and scored with the *true* parameters.  Regret is the true failure of the
    chosen scheme minus the true failure of the better one.
    """
    rng = np.random.default_rng(seed)
    regret, wrong = [], []
    for p2 in p2_cells:
        for tau in tau_cells:
            true = base.with_(p2=float(p2), tau_m=float(tau))
            t_mb, t_mf = measurement_based_fail(true), measurement_free_fail(true)
            f = np.exp(sigma * rng.standard_normal(4))
            believed = true.with_(p2=min(true.p2 * f[0], 0.3), eps0=min(true.eps0 * f[1], 0.2),
                                  p_meas=min(true.p_meas * f[2], 0.3), tau_m=true.tau_m * f[3])
            pick_mf = measurement_free_fail(believed) < measurement_based_fail(believed)
            chosen, best = (t_mf if pick_mf else t_mb), min(t_mb, t_mf)
            regret.append(chosen - best)
            wrong.append(float((t_mf < t_mb) != pick_mf))
    return {"regret": np.array(regret), "wrong": np.array(wrong)}


def params_dict(prm: Params) -> dict[str, float]:
    """Plain dict of parameters for result files."""
    return {k: (float(v) if not isinstance(v, bool) else v) for k, v in asdict(prm).items()}
