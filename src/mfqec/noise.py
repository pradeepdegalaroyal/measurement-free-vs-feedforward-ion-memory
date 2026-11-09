"""Exact bit-flip probability-tensor primitives.

A register of ``n`` bits is stored as a real tensor ``P`` of shape ``(2,) * n`` whose
entry ``P[b0, ..., b_{n-1}]`` is the probability of the classical bit string.  Because
the logical readout is a computational-basis measurement and the circuits contain only
CNOT / Toffoli permutations, a Pauli ``Z`` fault never turns into a bit flip, so the
``X``-component of every Pauli channel is a *complete* description of the readout
statistics.  All channels below therefore act as exact stochastic maps on ``P``.
(This holds for permutation gates treated as single units; see ``docs/METHOD.md`` for
the check against an explicit Hadamard/T decomposition of the Toffoli gate.)
"""
from __future__ import annotations

import numpy as np

Tensor = np.ndarray


def bit_flip(p: Tensor, axis: int, prob: float) -> Tensor:
    """Independent bit-flip channel with probability ``prob`` on one bit."""
    if prob <= 0.0:
        return p
    return (1.0 - prob) * p + prob * np.flip(p, axis=axis)


def depolarising2(p: Tensor, a: int, b: int, prob: float) -> Tensor:
    """Two-qubit depolarising channel of strength ``prob`` (X-component marginal).

    Of the 15 non-identity Paulis on qubits ``(a, b)``, 3 have no X or Y part (no flip),
    and 4 flip each of the patterns ``(a)``, ``(b)`` and ``(a, b)``.
    """
    if prob <= 0.0:
        return p
    fa = np.flip(p, axis=a)
    mixed = (3.0 * p + 4.0 * fa + 4.0 * np.flip(p, axis=b) + 4.0 * np.flip(fa, axis=b)) / 15.0
    return (1.0 - prob) * p + prob * mixed


def controlled_flip(p: Tensor, controls: tuple[tuple[int, int], ...], target: int) -> Tensor:
    """Flip ``target`` where every ``(axis, value)`` control is satisfied.

    ``controls=((c, 1),)`` is a CNOT, two controls give a Toffoli, and a control value
    of 0 gives a negatively-controlled gate (realised by free single-qubit X frames).
    """
    idx: list[object] = [slice(None)] * p.ndim
    for ax, val in controls:
        idx[ax] = val
    sub = p[tuple(idx)]
    shift = sum(1 for ax, _ in controls if ax < target)
    out = p.copy()
    out[tuple(idx)] = np.flip(sub, axis=target - shift)
    return out


def idle_prob(eps0: float, duration: float) -> float:
    """Bit-flip probability after ``duration`` gate-times of dephasing-type idling.

    ``eps0`` is the flip probability per gate-time; the flip probability follows
    ``(1 - exp(-2 gamma t)) / 2`` with ``gamma`` fixed by ``eps0``.
    """
    if eps0 <= 0.0 or duration <= 0.0:
        return 0.0
    gamma = -0.5 * np.log(1.0 - 2.0 * eps0)
    return float(0.5 * (1.0 - np.exp(-2.0 * gamma * duration)))


def combine(e1: float, e2: float) -> float:
    """Probability that exactly one of two independent flips happens."""
    return e1 * (1.0 - e2) + e2 * (1.0 - e1)


def append_zero_bits(p: Tensor, count: int) -> Tensor:
    """Append ``count`` fresh bits in state 0 as new trailing axes."""
    for _ in range(count):
        p = np.stack([p, np.zeros_like(p)], axis=-1)
    return p
