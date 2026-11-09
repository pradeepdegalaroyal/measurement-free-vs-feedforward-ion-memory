# Research gap

## State of the art (from the supplied literature)

The brief groups all thirty Scopus records under a single role, so the roles below are assigned by the author from publication date, citation count and article type.

* **Feed-forward quantum error correction (QEC) is the accepted baseline.** Surface-code and related memories measure syndromes mid-circuit, reset the ancillas and apply classically controlled corrections (Acharya et al., 2023, digital object identifier (DOI) 10.1038/s41586-022-05434-1; Zobrist et al., 2025, doi:10.1038/s41586-024-08449-y; Sivak et al., 2023, doi:10.1038/s41586-023-05782-6; Bluvstein et al., 2024, doi:10.1038/s41586-023-06927-3). On trapped ions, adaptive circuits with feed-forward have produced non-Abelian topological order on a 27-qubit state (Iqbal et al., 2024, doi:10.1038/s41586-023-06934-4).
* **Mid-circuit measurement is the stated pain point on ions.** The most recent trapped-ion paper in the set argues that feed-forward measurements are slow, hard to control and comparatively error-prone, and demonstrates universal logical operations on error-detecting codes with no mid-circuit measurement at all (Butt et al., 2026, doi:10.1038/s41467-026-68533-x).
* **Reviews point to measurement cost and error budgets.** The error-mitigation review lists open problems in reaching useful devices (Cai et al., 2023, doi:10.1103/RevModPhys.95.045005); the distributed-computing survey and the variational-algorithm review both ask how architectures behave as devices move towards fault tolerance (Caleffi et al., 2024, doi:10.1016/j.comnet.2024.110672; Wang et al., 2027, doi:10.15302/frontphys.2027.023301). A 2027 paper observes that measurements are costly and error-prone and that the intuition of omitting irrelevant ones has not been formalised (Matsumoto, 2027, doi:10.1007/s40314-026-03881-4).

## Related work table

| Paper | Year | Role (author-assigned) | What it does | What it does NOT cover |
|---|---|---|---|---|
| Butt et al., doi:10.1038/s41467-026-68533-x | 2026 | recent | Measurement-free universal logical gates on trapped-ion error-detecting codes | A quantitative break-even against feed-forward correction as a function of measurement duration and gate error |
| Matsumoto, doi:10.1007/s40314-026-03881-4 | 2027 | recent | Formalises which stabilisers are needed to correct erasures | Timing and error cost of the measurements; coherent alternatives |
| Singh et al., doi:10.1038/s41534-025-01146-2 | 2026 | recent | Thresholds of a distributed surface code versus entanglement scheme | Measurement-free correction; single-module ion memories |
| Joo, doi:10.1140/epjqt/s40507-026-00492-0 | 2026 | recent | Transformer decoder for the surface code | Hardware timing, ion-specific idle noise |
| Colmenarez et al., doi:10.1038/s41534-026-01362-4 | 2026 | recent | Unitary, distance-preserving logical state preparation | Repeated correction cycles |
| Acharya et al., doi:10.1038/s41586-022-05434-1 | 2023 | landmark | Surface-code logical-qubit performance across code sizes | Trapped-ion timing; measurement-free variants |
| Zobrist et al., doi:10.1038/s41586-024-08449-y | 2025 | landmark | Below-threshold surface-code memories with a real-time decoder | Cost of slow measurement; coherent correction |
| Sivak et al., doi:10.1038/s41586-023-05782-6 | 2023 | landmark | Real-time QEC beyond break-even (feed-forward baseline) | Ion platforms; any measurement-free comparison |
| Bluvstein et al., doi:10.1038/s41586-023-06927-3 | 2024 | landmark | Logical processor on neutral atoms with mid-circuit readout | Cost trade-off of removing mid-circuit readout |
| Iqbal et al., doi:10.1038/s41586-023-06934-4 | 2024 | landmark | Adaptive trapped-ion circuits with feed-forward | Error-correction memory comparison |
| Cai et al., doi:10.1103/RevModPhys.95.045005 | 2023 | review | Survey of error mitigation and open problems | Measurement-free correction budgets |
| Caleffi et al., doi:10.1016/j.comnet.2024.110672 | 2024 | review | Survey of distributed quantum computing | Single-module correction cycle design |
| Bartolucci et al., doi:10.1038/s41467-023-36493-1 | 2023 | adjacent | Fusion-based computation, measurement-centred fault tolerance | Matter-qubit memories with idle noise |
| van den Berg et al., doi:10.1038/s41567-023-02042-2 | 2023 | adjacent | Probabilistic error cancellation with sparse noise models | Encoded memories; scheme selection from noise estimates |

## The gap

Feed-forward correction is the baseline in every landmark experiment above, and the newest trapped-ion work removes mid-circuit measurement entirely. What is missing is a **quantitative, noise-parametrised answer to a design question an ion-trap experimenter faces before booking machine time: for which (measurement duration, two-qubit error, idle dephasing, readout error, Toffoli cost) is coherent, measurement-free correction actually better than measurement plus feed-forward, and how many shots does it take to find out?** Measurement-free demonstrations report that the approach is feasible; they do not map where it wins.

**Novelty (one sentence).** This study supplies exact break-even maps between measurement-free and feed-forward correction of a repetition memory on a trapped-ion-style noise model, together with the shot budget and noise-characterisation tolerance needed to choose between them.

## Hypotheses and success metrics

Scale: exact classical-probability simulation of a three-ion bit-flip repetition memory with two ancilla ions, one to five cycles (up to 13 ions), on a laptop. Parameters are illustrative and not calibrated to IonQ Aria or Forte.

* **H1 (window exists).** With a decomposed Toffoli (six two-qubit gates), idle flip 1e-3 and readout flip 1e-3 per gate-time, measurement-free correction beats feed-forward above a break-even measurement time tau*(p2) that is at most 100 two-qubit gate-times for every two-qubit error p2 in [1e-4, 1e-2]. Falsified if tau* > 100 anywhere in that range.
* **H2 (single-ratio collapse).** The excess tau* - 3 n2 depends on p2 / eps0 only: log-log slope in [0.8, 1.2] and scatter of at most 0.15 decades across five idle rates eps0.
* **H3 (readout is not enough).** Removing mid-circuit readout error alone does not make measurement-free correction win for instantaneous measurement (tau_m = 2) unless the readout flip exceeds 1e-2 for the exact Toffoli.
* **H4 (shot cost).** The correct scheme can be identified with at most 20,000 shots per scheme for every tau_m outside a band of +/-10 % around the break-even.
* **H5 (calibration tolerance).** Choosing the scheme from a noise model whose parameters carry a log-normal error of spread 0.3 mis-classifies at most 10 % of the (p2, tau_m) plane.

H1 and H2 thresholds were fixed after a coarse pilot scan, so they are exploratory; H3-H5 were fixed before the corresponding experiments.
