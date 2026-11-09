"""Three-bit repetition memory: measurement-based vs measurement-free correction.

Both schemes protect one logical bit stored in three data ions against dephasing-type
errors (bit flips in the Hadamard frame).  Each cycle extracts the two parities
``d0+d1`` and ``d1+d2`` into two ancilla ions and then corrects the data:

* **feed-forward** (``measurement_based``): ancillas are measured (cost ``tau_m``
  gate-times, readout flip ``p_meas``, spectator error ``p_xt``), reset, and a classically
  controlled X is applied.
* **measurement-free** (``measurement_free``): fresh ancillas every cycle, three
  coherently controlled flips (each Toffoli costs ``n2`` noisy two-qubit gates); no
  mid-circuit measurement and no reset.

Time is measured in units of one two-qubit gate duration.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .noise import (Tensor, append_zero_bits, bit_flip, combine, controlled_flip,
                    depolarising2, idle_prob)

SYNDROME_GATES = ((0, 3), (1, 3), (1, 4), (2, 4))  # (data, ancilla) CNOT pairs
SYNDROME_TIME = float(len(SYNDROME_GATES))
_CORRECTION = {(1, 1): 1, (1, 0): 0, (0, 1): 2}     # reported syndrome -> data bit to flip


@dataclass(frozen=True)
class Params:
    """Physical noise and timing parameters (all times in two-qubit gate-times)."""

    p2: float = 5e-3          # two-qubit depolarising strength
    eps0: float = 1e-3        # idle bit-flip probability per gate-time
    p_meas: float = 5e-3      # readout flip probability (mid-circuit and final)
    tau_m: float = 10.0       # mid-circuit measurement + reset + feed-forward time
    p_xt: float = 0.0         # spectator flip on each data ion per mid-circuit measurement
    n2: int = 6               # two-qubit gates per coherent Toffoli
    p1_ratio: float = 0.1     # single-qubit gate error as a fraction of p2
    rounds: int = 3           # correction cycles before the final majority readout

    def with_(self, **kw: float) -> "Params":
        """Return a copy with some fields replaced."""
        return replace(self, **kw)


def _syndrome(p: Tensor, prm: Params, a1: int = 3) -> Tensor:
    """Four syndrome CNOTs onto the ancilla pair stored at axes ``(a1, a1 + 1)``."""
    for d, a in SYNDROME_GATES:
        anc = a1 + (a - 3)
        p = controlled_flip(p, ((d, 1),), anc)
        p = depolarising2(p, d, anc, prm.p2)
    return p


def _idle_data(p: Tensor, prm: Params, duration: float, extra: float = 0.0) -> Tensor:
    e = combine(idle_prob(prm.eps0, duration), extra)
    for q in range(3):
        p = bit_flip(p, q, e)
    return p


def _marginal_data(p: Tensor) -> Tensor:
    return p.sum(axis=tuple(range(3, p.ndim))) if p.ndim > 3 else p


def _readout_fail(p: Tensor, prm: Params) -> float:
    """Probability that the noisy majority vote of the three data bits is wrong."""
    d = _marginal_data(p)
    for q in range(3):
        d = bit_flip(d, q, prm.p_meas)
    ones = np.add.outer(np.add.outer(np.arange(2), np.arange(2)), np.arange(2))
    return float(d[ones >= 2].sum())


def measurement_based_fail(prm: Params) -> float:
    """Logical failure probability of the feed-forward memory (exact)."""
    data = np.zeros((2, 2, 2))
    data[0, 0, 0] = 1.0
    flip1 = (2.0 / 3.0) * prm.p1_ratio * prm.p2
    for _ in range(prm.rounds):
        p = append_zero_bits(data, 2)
        p = _syndrome(p, prm)
        p = _idle_data(p, prm, SYNDROME_TIME)
        p = _idle_data(p, prm, prm.tau_m, prm.p_xt)  # measurement window
        branches = {s: np.zeros((2, 2, 2)) for s in ((0, 0), (0, 1), (1, 0), (1, 1))}
        for a1 in (0, 1):
            for a2 in (0, 1):
                d_a = p[:, :, :, a1, a2]
                for m1 in (0, 1):
                    for m2 in (0, 1):
                        w = ((prm.p_meas if m1 else 1 - prm.p_meas)
                             * (prm.p_meas if m2 else 1 - prm.p_meas))
                        branches[(a1 ^ m1, a2 ^ m2)] += w * d_a
        data = np.zeros((2, 2, 2))
        for s, d_s in branches.items():
            if s in _CORRECTION:
                d_s = (1 - flip1) * np.flip(d_s, axis=_CORRECTION[s]) + flip1 * d_s
            data += d_s
    return _readout_fail(data, prm)


def measurement_free_fail(prm: Params) -> float:
    """Logical failure probability of the coherent-correction memory (exact)."""
    p = np.zeros((2, 2, 2))
    p[0, 0, 0] = 1.0
    for _ in range(prm.rounds):
        a1, a2 = p.ndim, p.ndim + 1
        p = append_zero_bits(p, 2)
        p = _syndrome(p, prm, a1)
        p = _idle_data(p, prm, SYNDROME_TIME)
        for v1, v2, tgt in ((1, 1, 1), (1, 0, 0), (0, 1, 2)):
            p = controlled_flip(p, ((a1, v1), (a2, v2)), tgt)
            pairs = ((a1, tgt), (a2, tgt), (a1, a2))
            for k in range(prm.n2):
                x, y = pairs[k % 3]
                p = depolarising2(p, x, y, prm.p2)
            p = _idle_data(p, prm, float(prm.n2))
    return _readout_fail(p, prm)


def total_time_mb(prm: Params) -> float:
    """Wall-clock duration (gate-times) of the feed-forward memory."""
    return prm.rounds * (SYNDROME_TIME + prm.tau_m)


def total_time_mf(prm: Params) -> float:
    """Wall-clock duration (gate-times) of the measurement-free memory."""
    return prm.rounds * (SYNDROME_TIME + 3 * prm.n2)


def bare_fail(prm: Params, duration: float) -> float:
    """Unencoded qubit idling for ``duration`` gate-times then read out once."""
    return combine(idle_prob(prm.eps0, duration), prm.p_meas)


def qubit_count(prm: Params, scheme: str) -> int:
    """Number of physical ions: fresh ancillas every cycle for measurement-free."""
    return 3 + (2 * prm.rounds if scheme == "mf" else 2)
