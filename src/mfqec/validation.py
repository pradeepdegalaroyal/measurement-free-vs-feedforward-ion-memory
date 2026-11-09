"""Cross-checks of the classical bit-flip model against independent simulators.

* :func:`aer_measurement_free_fail` builds the one-cycle measurement-free circuit with an
  explicit Hadamard/T/CNOT Toffoli decomposition and simulates it with Qiskit Aer
  density matrices.  This tests the one assumption of the fast model: that a Toffoli can
  be treated as an ideal permutation followed by ``n2`` pair-depolarising channels.
* :func:`trajectory_fail` is an independent bit-level Monte Carlo sampler for both schemes.
"""
from __future__ import annotations

import numpy as np

from .schemes import SYNDROME_GATES, SYNDROME_TIME, Params


def aer_measurement_free_fail(p2: float, decompose: bool = True) -> float:
    """Exact one-cycle logical failure of the measurement-free memory in Qiskit Aer.

    Only two-qubit depolarising noise on every CNOT (strength ``p2``); no idle or readout
    noise.  The failure is the probability that at least two of the three data bits read 1.
    """
    from qiskit import QuantumCircuit, transpile
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error

    qc = QuantumCircuit(5)
    for d, a in SYNDROME_GATES:
        qc.cx(d, a)
    for v1, v2, tgt in ((1, 1, 1), (1, 0, 0), (0, 1, 2)):
        if v1 == 0:
            qc.x(3)
        if v2 == 0:
            qc.x(4)
        qc.ccx(3, 4, tgt)
        if v1 == 0:
            qc.x(3)
        if v2 == 0:
            qc.x(4)
    qc.save_probabilities(qubits=[0, 1, 2], label="pd")
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(depolarizing_error(p2 * 16.0 / 15.0, 2), ["cx"])
    sim = AerSimulator(method="density_matrix", noise_model=nm)
    basis = ["cx", "h", "t", "tdg", "x", "save_probabilities"]
    tqc = transpile(qc, sim, basis_gates=basis, optimization_level=0)
    probs = sim.run(tqc).result().data()["pd"]
    arr = np.asarray(probs, dtype=float)
    return float(sum(arr[i] for i in range(8) if bin(i).count("1") >= 2))


def trajectory_fail(prm: Params, scheme: str, shots: int, seed: int) -> float:
    """Independent bit-level Monte Carlo estimate of the logical failure probability."""
    rng = np.random.default_rng(seed)
    gamma = -0.5 * np.log(1 - 2 * prm.eps0) if prm.eps0 > 0 else 0.0

    def idle_p(t: float) -> float:
        return 0.5 * (1 - np.exp(-2 * gamma * t))

    def comb(a: float, b: float) -> float:
        return a * (1 - b) + b * (1 - a)

    def dep(bits: np.ndarray, a: int, b: int) -> None:
        hit = rng.random(shots) < prm.p2
        pick = rng.integers(1, 16, size=shots)  # 4-bit Pauli index (xa, za, xb, zb)
        xa = (pick >> 3) & 1
        xb = (pick >> 1) & 1
        bits[:, a] ^= (hit & (xa == 1)).astype(np.uint8)
        bits[:, b] ^= (hit & (xb == 1)).astype(np.uint8)

    def idle(bits: np.ndarray, e: float) -> None:
        for q in range(3):
            bits[:, q] ^= (rng.random(shots) < e).astype(np.uint8)

    data = np.zeros((shots, 3), dtype=np.uint8)
    flip1 = (2.0 / 3.0) * prm.p1_ratio * prm.p2
    for _ in range(prm.rounds):
        bits = np.concatenate([data, np.zeros((shots, 2), dtype=np.uint8)], axis=1)
        for d, a in SYNDROME_GATES:
            bits[:, a] ^= bits[:, d]
            dep(bits, d, a)
        idle(bits, idle_p(SYNDROME_TIME))
        if scheme == "mb":
            idle(bits, comb(idle_p(prm.tau_m), prm.p_xt))
            s1 = bits[:, 3] ^ (rng.random(shots) < prm.p_meas)
            s2 = bits[:, 4] ^ (rng.random(shots) < prm.p_meas)
            for v1, v2, tgt in ((1, 1, 1), (1, 0, 0), (0, 1, 2)):
                m = (s1 == v1) & (s2 == v2)
                err = rng.random(shots) < flip1
                bits[:, tgt] ^= (m & ~err).astype(np.uint8)
            bits[:, 3:] = 0
            data = bits[:, :3]
        else:
            a1, a2 = 3, 4
            for v1, v2, tgt in ((1, 1, 1), (1, 0, 0), (0, 1, 2)):
                m = (bits[:, a1] == v1) & (bits[:, a2] == v2)
                bits[:, tgt] ^= m.astype(np.uint8)
                pairs = ((a1, tgt), (a2, tgt), (a1, a2))
                for k in range(prm.n2):
                    dep(bits, *pairs[k % 3])
                idle(bits, idle_p(float(prm.n2)))
            data = bits[:, :3]
    for q in range(3):
        data[:, q] ^= (rng.random(shots) < prm.p_meas).astype(np.uint8)
    return float((data.sum(axis=1) >= 2).mean())
