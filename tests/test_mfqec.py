"""Fast tests for the exact bit-flip model."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mfqec import (Params, bare_fail, measurement_based_fail, measurement_free_fail,  # noqa: E402
                   qubit_count, total_time_mb, total_time_mf)
from mfqec.analysis import shots_to_resolve, tau_star  # noqa: E402
from mfqec.noise import controlled_flip, depolarising2, idle_prob  # noqa: E402
from mfqec.validation import trajectory_fail  # noqa: E402


def test_noiseless_memory_never_fails():
    prm = Params(p2=0.0, eps0=0.0, p_meas=0.0, p_xt=0.0)
    assert measurement_based_fail(prm) == 0.0
    assert measurement_free_fail(prm) == 0.0


def test_channels_conserve_probability():
    rng = np.random.default_rng(0)
    p = rng.random((2,) * 5)
    p /= p.sum()
    q = depolarising2(p, 1, 4, 0.2)
    assert q.sum() == pytest.approx(1.0)
    r = controlled_flip(q, ((0, 1), (3, 0)), 2)
    assert r.sum() == pytest.approx(1.0)


def test_toffoli_permutation_is_correct():
    p = np.zeros((2,) * 3)
    p[1, 1, 0] = 1.0
    out = controlled_flip(p, ((0, 1), (1, 1)), 2)
    assert out[1, 1, 1] == 1.0
    p2 = np.zeros((2,) * 3)
    p2[1, 0, 0] = 1.0
    assert controlled_flip(p2, ((0, 1), (1, 1)), 2)[1, 0, 0] == 1.0


def test_idle_probability_matches_definition():
    assert idle_prob(1e-3, 1.0) == pytest.approx(1e-3)
    assert idle_prob(1e-3, 0.0) == 0.0
    assert idle_prob(1e-3, 1e6) == pytest.approx(0.5, abs=1e-6)


def test_feedforward_gets_worse_with_measurement_time_and_mf_does_not_depend_on_it():
    base = Params()
    fast, slow = base.with_(tau_m=1.0), base.with_(tau_m=100.0)
    assert measurement_based_fail(slow) > measurement_based_fail(fast)
    assert measurement_free_fail(slow) == measurement_free_fail(fast)


def test_bare_qubit_with_only_readout_noise():
    prm = Params(eps0=0.0, p_meas=0.01)
    assert bare_fail(prm, 50.0) == pytest.approx(0.01)


def test_break_even_is_a_crossing_and_tends_to_idle_balance():
    prm = Params(p2=1e-4, eps0=1e-3, p_meas=1e-4, n2=6, rounds=3)
    ts = tau_star(prm)
    assert measurement_based_fail(prm.with_(tau_m=ts * 0.8)) < measurement_free_fail(prm)
    assert measurement_based_fail(prm.with_(tau_m=ts * 1.25)) > measurement_free_fail(prm)
    assert 15.0 < ts < 25.0  # near 3 * n2 = 18 gate-times when gates are nearly perfect


@pytest.mark.parametrize("scheme,fn", [("mb", measurement_based_fail), ("mf", measurement_free_fail)])
def test_exact_model_matches_independent_trajectories(scheme, fn):
    prm = Params(p2=1e-2, eps0=2e-3, p_meas=5e-3, tau_m=8.0, p_xt=2e-3, rounds=2)
    shots = 300_000
    est = trajectory_fail(prm, scheme, shots, 7)
    exact = fn(prm)
    sigma = np.sqrt(exact * (1 - exact) / shots)
    assert abs(est - exact) < 5 * sigma


def test_cost_bookkeeping():
    prm = Params(rounds=3, n2=6, tau_m=10.0)
    assert qubit_count(prm, "mb") == 5 and qubit_count(prm, "mf") == 9
    assert total_time_mb(prm) == 3 * 14 and total_time_mf(prm) == 3 * 22


def test_shots_to_resolve_scaling():
    assert shots_to_resolve(0.01, 0.01) == float("inf")
    assert shots_to_resolve(0.01, 0.02) < shots_to_resolve(0.01, 0.012)
