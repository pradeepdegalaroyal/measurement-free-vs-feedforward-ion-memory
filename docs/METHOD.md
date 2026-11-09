# Method

## System

One logical bit is stored in three data ions (bit-flip repetition code; in the Hadamard frame this is a dephasing-protected memory, which is the dominant idle error of ion qubits). Two ancilla ions extract the parities d0+d1 and d1+d2 with four controlled-NOT (CNOT) gates per cycle. After `R` cycles the three data ions are read out and the logical bit is the majority vote. Time is counted in units of one two-qubit gate duration.

* **Feed-forward (standard quantum baseline from the literature).** Ancillas are measured (duration `tau_m`, readout flip `p_meas`, optional spectator flip `p_xt` on each data ion), reset, and a lookup table applies a classically controlled X to the implicated data ion. Two reusable ancillas.
* **Measurement-free (the alternative under test).** Fresh ancillas each cycle; three coherent Toffoli gates (ancilla pattern 11 flips d1, 10 flips d0, 01 flips d2; negative controls come from free single-qubit X frames). Each Toffoli costs `n2` two-qubit gates (6 exact, 3 relative-phase, 1 native). No measurement, no reset, `2R + 3` ions.
* **Classical baseline.** An unencoded qubit idling for the same wall-clock duration and read out once. A lookup-table classical decoder is part of the feed-forward arm.

## Noise model

Two-qubit depolarising error `p2` after every two-qubit gate (Toffoli: `n2` such channels on the three qubit pairs in turn), single-qubit gate error `0.1 p2`, idle bit-flip `eps0` per gate-time following (1 - exp(-2 gamma t)) / 2, readout flip `p_meas` on every measurement including the final readout, and `p_xt` spectator flips during mid-circuit measurement (feed-forward only).

## Exact simulation

Because the readout is in the computational basis and the circuits contain only permutation gates, a Z fault never becomes a bit flip; the X-component of each Pauli channel therefore gives the exact output distribution. The register is a probability tensor over bit strings, and every channel is an exact stochastic map. Logical failure probabilities are exact (no sampling) and cost milliseconds, enabling dense maps.

Checks: (i) one-cycle measurement-free circuit with an explicit Hadamard/T/CNOT Toffoli decomposition simulated with Qiskit Aer density matrices agrees with the model; (ii) an independent bit-level Monte Carlo sampler agrees for both schemes over three cycles; (iii) unit tests on probability conservation and limits.

## Metrics

* Logical failure probability after `R` cycles; break-even measurement time `tau*` (bisection on log `tau_m`; feed-forward failure is monotone in `tau_m`, measurement-free does not depend on it); log10 ratio of the two failure probabilities for the phase maps.
* Pseudo-threshold: largest `p2` at which a scheme beats the unencoded qubit of equal duration.
* Shots per scheme to resolve the winner (two-proportion test, 95 % level, 80 % power), and a five-seed binomial emulation of a 20,000-shot experiment.

## Ablation, robustness, seeds

* **Ablation:** gates per Toffoli `n2` in {1, 2, 3, 6, 9}; number of cycles 1-5.
* **Noise and robustness sweeps:** idle flip `eps0` 1e-4 to 1e-2, readout flip 1e-4 to 3e-2, spectator flip 0 to 1e-2, readout threshold at tau_m = 2, and a noise-characterisation-error study (log-normal spread 0, 0.15, 0.3, 0.6 on `p2`, `eps0`, `p_meas`, `tau_m`).
* **Seeds:** the model is exact; seeds (0-4) govern finite-shot emulation, trajectory cross-checks and characterisation-error draws. Statistics are mean +/- sample standard deviation over five seeds.

## Known limitation of the model

The Toffoli is treated as an ideal permutation followed by `n2` pair-depolarising channels; the Aer check confirms this for one cycle at the tested `p2`. Only the Z-basis logical value is scored, so phase coherence of the logical qubit is not assessed.
