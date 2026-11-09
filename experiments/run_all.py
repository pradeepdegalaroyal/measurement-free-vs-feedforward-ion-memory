"""Reproduce every result and figure.  ``python experiments/run_all.py [--quick]``."""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt  # noqa: E402

from mfqec import (Params, bare_fail, measurement_based_fail, measurement_free_fail,  # noqa: E402
                   qubit_count, total_time_mb, total_time_mf)
from mfqec.analysis import (calibration_regret, finite_shot_difference, phase_map,  # noqa: E402
                            shots_to_resolve, tau_star)
from mfqec.plots import colors, phase_panels, save  # noqa: E402
from mfqec.validation import aer_measurement_free_fail, trajectory_fail  # noqa: E402

BASE = Params(p2=1e-3, eps0=1e-3, p_meas=1e-3, tau_m=10.0, n2=6, rounds=3)
SEEDS = [0, 1, 2, 3, 4]


def timed(name: str, t0: float) -> None:
    print(f"[{time.time() - t0:6.1f}s] {name}", flush=True)


def exp_validation(res: Path, quick: bool) -> dict:
    """Aer density-matrix and bit-level trajectory cross-checks of the fast model."""
    rows = []
    for p2 in (1e-3, 5e-3, 2e-2):
        prm = Params(p2=p2, eps0=0.0, p_meas=0.0, rounds=1)
        rows.append({"check": "aer_explicit_toffoli_1cycle", "p2": p2,
                     "model": measurement_free_fail(prm), "reference": aer_measurement_free_fail(p2)})
    shots = 100_000 if quick else 1_000_000
    prm = Params(p2=3e-3, eps0=1e-3, p_meas=2e-3, tau_m=15.0, p_xt=1e-3, rounds=3)
    for scheme, fn in (("mb", measurement_based_fail), ("mf", measurement_free_fail)):
        vals = [trajectory_fail(prm, scheme, shots, 10 + s) for s in SEEDS]
        m = fn(prm)
        rows.append({"check": f"trajectory_{scheme}_3cycles", "p2": prm.p2, "model": m,
                     "reference": float(np.mean(vals)), "reference_std": float(np.std(vals, ddof=1)),
                     "shots_per_seed": shots})
    df = pd.DataFrame(rows)
    df.to_csv(res / "validation.csv", index=False)
    return {"max_abs_diff_aer": float((df[df.check.str.startswith("aer")].eval("abs(model-reference)")).max())}


def exp_phase_maps(res: Path, fig: Path, quick: bool) -> dict:
    """E1: phase diagrams over (p2, tau_m) for three Toffoli costs."""
    npt = 8 if quick else 36
    p2g = np.logspace(-4, -2, npt)
    tg = np.logspace(np.log10(0.5), np.log10(200), npt)
    rows, maps = [], {}
    for n2 in (6, 3, 1):
        out = phase_map(BASE.with_(n2=n2), p2g, tg)
        maps[n2] = np.log10(out["mb"] / out["mf"])
        for i, t in enumerate(tg):
            for j, p in enumerate(p2g):
                rows.append({"n2": n2, "p2": p, "tau_m": t, "mb": out["mb"][i, j],
                             "mf": out["mf"][i, j], "log10_ratio": maps[n2][i, j]})
    pd.DataFrame(rows).to_csv(res / "phase_maps.csv", index=False)
    phase_panels(maps, p2g, tg, fig / "fig1_phase_maps.png")
    return {f"mf_win_fraction_n2_{n2}": float((m > 0).mean()) for n2, m in maps.items()}


def exp_tau_star(res: Path, fig: Path, quick: bool) -> dict:
    """E2: break-even measurement time vs p2 for ablation (n2), idle noise, readout noise."""
    rows = []
    P2_AXIS = np.logspace(-4, -2, 5 if quick else 13)
    for n2 in (1, 2, 3, 6, 9):
        for p2 in P2_AXIS:
            rows.append({"sweep": "toffoli_cost", "value": n2, "p2": p2,
                         "tau_star": tau_star(BASE.with_(p2=p2, n2=n2))})
    for e in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
        for p2 in P2_AXIS:
            rows.append({"sweep": "idle_eps0", "value": e, "p2": p2,
                         "tau_star": tau_star(BASE.with_(p2=p2, eps0=e))})
    for pm in (1e-4, 1e-3, 1e-2, 3e-2):
        for p2 in P2_AXIS:
            rows.append({"sweep": "readout_p_meas", "value": pm, "p2": p2,
                         "tau_star": tau_star(BASE.with_(p2=p2, p_meas=pm))})
    for px in (0.0, 1e-4, 1e-3, 1e-2):
        for p2 in P2_AXIS:
            rows.append({"sweep": "spectator_p_xt", "value": px, "p2": p2,
                         "tau_star": tau_star(BASE.with_(p2=p2, p_xt=px))})
    df = pd.DataFrame(rows)
    df.to_csv(res / "tau_star.csv", index=False)

    fig_, axes = plt.subplots(1, 4, figsize=(17, 4.0), sharey=True)
    titles = {"toffoli_cost": "gates per Toffoli $n_2$", "idle_eps0": "idle flip per gate-time $\\varepsilon_0$",
              "readout_p_meas": "readout flip $p_{\\mathrm{meas}}$", "spectator_p_xt": "spectator flip $p_{xt}$"}
    mk = ["o", "s", "^", "D", "v"]
    ls = ["-", "--", "-.", ":", (0, (5, 1))]
    for ax, (name, g) in zip(axes, df.groupby("sweep", sort=False)):
        vals = list(dict.fromkeys(g.value))
        for k, (v, c) in enumerate(zip(vals, colors(len(vals)))):
            s = g[g.value == v]
            ax.plot(s.p2, s.tau_star, color=c, marker=mk[k % 5], ls=ls[k % 5], ms=4.5, lw=1.8,
                    label=f"{v:g}")
        ax.axhline(0.5, color="0.6", lw=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("two-qubit error probability $p_2$")
        ax.set_title(titles[name], fontsize=10)
        ax.legend(fontsize=8, title=None)
    axes[0].set_ylabel("break-even $\\tau_m^*$ (gate-times)")
    fig_.suptitle("Break-even measurement time: measurement-free wins above each curve", y=1.0)
    save(fig_, fig / "fig2_break_even_sweeps.png")

    # H2: does tau* - 3 n2 collapse on p2 / eps0 ?
    e = df[df.sweep == "idle_eps0"].copy()
    e["excess"] = e.tau_star - 18.0
    e["x"] = e.p2 / e.value
    e = e[(e.excess > 0.5) & (e.tau_star < 9_000)]
    slope = float(np.polyfit(np.log(e.x), np.log(e.excess), 1)[0])
    spread = float(np.std(np.log10(e.excess) - np.polyval(np.polyfit(np.log10(e.x), np.log10(e.excess), 1),
                                                          np.log10(e.x))))
    tb = df[(df.sweep == "toffoli_cost")]
    low = tb[tb.p2 == tb.p2.min()].set_index("value").tau_star
    return {"collapse_loglog_slope": slope, "collapse_resid_dex": spread,
            "tau_star_at_min_p2": {str(int(k)): float(v) for k, v in low.items()}}


def exp_finite_shots(res: Path, fig: Path, quick: bool) -> dict:
    """E3: bare vs both schemes vs tau_m, with a finite-shot (hardware-like) emulation."""
    shots = 20_000
    taus = np.logspace(np.log10(0.5), np.log10(150), 8 if quick else 16)
    rows = []
    mf = measurement_free_fail(BASE)
    for t in taus:
        prm = BASE.with_(tau_m=float(t))
        mb = measurement_based_fail(prm)
        fs = finite_shot_difference(mb, mf, shots, SEEDS)
        rows.append({"tau_m": t, "mb": mb, "mf": mf, "bare_same_time_mb": bare_fail(prm, total_time_mb(prm)),
                     "bare_same_time_mf": bare_fail(prm, total_time_mf(prm)),
                     "diff_exact": mb - mf, "diff_shots_mean": fs["diff"].mean(),
                     "diff_shots_std": fs["diff"].std(ddof=1), "correct_sign_seeds": int(
                         ((fs["diff"] > 0) == (mb > mf)).sum()),
                     "shots_to_resolve": shots_to_resolve(mb, mf)})
    df = pd.DataFrame(rows)
    df.to_csv(res / "finite_shot.csv", index=False)

    fig_, (a, b) = plt.subplots(1, 2, figsize=(11, 4.2))
    c = colors(4)
    rng_seed = [np.random.default_rng(2000 + s) for s in SEEDS]
    est_mb = np.array([[r.binomial(shots, m) / shots for m in df.mb] for r in rng_seed])
    est_mf = np.array([[r.binomial(shots, m) / shots for m in df.mf] for r in rng_seed])
    for est, col, mk, lab, ex in ((est_mb, c[0], "o", "feed-forward", df.mb),
                                  (est_mf, c[2], "s", "measurement-free", df.mf)):
        a.fill_between(df.tau_m, est.mean(0) - est.std(0, ddof=1), est.mean(0) + est.std(0, ddof=1),
                       color=col, alpha=0.3)
        a.plot(df.tau_m, est.mean(0), color=col, marker=mk, ms=4, ls="-", label=f"{lab}: mean $\\pm$ std, 5 seeds")
        a.plot(df.tau_m, ex, color=col, ls=":", lw=1)
    a.plot(df.tau_m, df.bare_same_time_mb, color=c[3], ls="--", lw=1.6, label="unencoded qubit, same duration")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("measurement time $\\tau_m$ (two-qubit gate-times)")
    a.set_ylabel(f"logical failure probability ({shots:,} shots)")
    a.set_title("Finite-shot estimate, $p_2$ = 1e-3, 3 cycles (dotted: exact)", fontsize=10)
    a.legend(fontsize=8)
    b.plot(df.tau_m, df.shots_to_resolve, color=c[1], marker="D", ms=4, ls="-.", lw=1.8)
    b.set_xscale("log")
    b.set_yscale("log")
    b.axhline(shots, color="0.5", ls=":", lw=1)
    b.text(df.tau_m.min() * 1.05, shots * 1.15, "20,000-shot budget", fontsize=8, color="0.3")
    b.set_xlabel("measurement time $\\tau_m$ (two-qubit gate-times)")
    b.set_ylabel("shots per scheme to resolve winner")
    b.set_title("Cost of deciding (95 % level, 80 % power)", fontsize=10)
    save(fig_, fig / "fig3_finite_shot_cost.png")
    near = df.iloc[(df.diff_exact.abs()).argmin()]
    return {"tau_cross_grid_nearest": float(near.tau_m), "shots_to_resolve_at_cross": float(near.shots_to_resolve),
            "min_shots_to_resolve": float(df.shots_to_resolve.min()),
            "correct_sign_fraction_outside_10pct_of_cross": float(
                df[(df.tau_m < 0.9 * near.tau_m) | (df.tau_m > 1.1 * near.tau_m)].correct_sign_seeds.mean() / 5)}


def exp_rounds(res: Path, fig: Path, quick: bool) -> dict:
    """E4: scaling with the number of correction cycles (ions and failure probability)."""
    rows = []
    for r in ([1, 2, 3] if quick else [1, 2, 3, 4, 5]):
        for tau in (5.0, 20.0, 60.0):
            prm = BASE.with_(rounds=r, tau_m=tau)
            rows.append({"rounds": r, "tau_m": tau, "mb": measurement_based_fail(prm),
                         "mf": measurement_free_fail(prm),
                         "bare_mb_time": bare_fail(prm, total_time_mb(prm)),
                         "bare_mf_time": bare_fail(prm, total_time_mf(prm)),
                         "ions_mb": qubit_count(prm, "mb"), "ions_mf": qubit_count(prm, "mf")})
    df = pd.DataFrame(rows)
    df.to_csv(res / "rounds.csv", index=False)
    fig_, ax = plt.subplots(figsize=(6.4, 4.4))
    c = colors(3)
    ls_mb = {5.0: "-", 20.0: "--", 60.0: "-."}
    for k, tau in enumerate((5.0, 20.0, 60.0)):
        s = df[df.tau_m == tau]
        ax.plot(s.rounds, s.mb, color=c[k], marker="o", ls=ls_mb[tau], label=f"feed-forward, $\\tau_m$ = {tau:g}")
    s = df[df.tau_m == 20.0]
    ax.plot(s.rounds, s.mf, color="black", marker="s", ls=":", lw=2, label="measurement-free (any $\\tau_m$)")
    ax.set_yscale("log")
    ax.set_xlabel("correction cycles")
    ax.set_xticks(sorted(df.rounds.unique()))
    ax.set_ylabel("logical failure probability")
    ax.set_title("Cycle scaling at $p_2$ = 1e-3")
    ax.legend(fontsize=8)
    save(fig_, fig / "fig4_cycle_scaling.png")
    r3 = df[(df.rounds == 3) & (df.tau_m == 20.0)].iloc[0]
    return {"ions_mb_3cycles": int(r3.ions_mb), "ions_mf_3cycles": int(r3.ions_mf)}


def exp_thresholds(res: Path) -> dict:
    """E5: readout error above which measurement-free wins even at tau_m = 2."""
    rows = []
    for n2 in (1, 3, 6):
        for p2 in (1e-4, 1e-3, 1e-2):
            lo, hi, found = 1e-5, 0.3, np.nan
            prm = BASE.with_(n2=n2, p2=p2, tau_m=2.0)
            f = lambda pm: measurement_based_fail(prm.with_(p_meas=pm)) - measurement_free_fail(prm.with_(p_meas=pm))
            if f(lo) > 0:
                found = lo
            elif f(hi) > 0:
                for _ in range(40):
                    mid = np.sqrt(lo * hi)
                    lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
                found = np.sqrt(lo * hi)
            rows.append({"n2": n2, "p2": p2, "p_meas_star": found})
    df = pd.DataFrame(rows)
    df.to_csv(res / "readout_threshold.csv", index=False)
    return {f"n2={r.n2},p2={r.p2:g}": (None if np.isnan(r.p_meas_star) else float(r.p_meas_star))
            for r in df.itertuples()}


def exp_calibration(res: Path, fig: Path, quick: bool) -> dict:
    """E6: regret of choosing the scheme from a mis-characterised noise model."""
    k = 4 if quick else 8
    p2c = np.logspace(-4, -2, k)
    tc = np.logspace(0, np.log10(200), k)
    sigmas = [0.0, 0.15, 0.3, 0.6]
    rows = []
    for s in sigmas:
        for seed in SEEDS:
            o = calibration_regret(BASE, p2c, tc, s, seed)
            rows.append({"sigma": s, "seed": seed, "misclassified": o["wrong"].mean(),
                         "mean_regret": o["regret"].mean()})
    df = pd.DataFrame(rows)
    df.to_csv(res / "calibration_regret.csv", index=False)
    fig_, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    c = colors(len(sigmas))
    for ax, col, lab in ((axes[0], "misclassified", "fraction of cells with wrong scheme"),
                         (axes[1], "mean_regret", "mean excess logical failure (regret)")):
        data = [df[df.sigma == s][col].values for s in sigmas]
        parts = ax.violinplot(data, positions=range(len(sigmas)), showextrema=False, widths=0.8)
        for body, cc in zip(parts["bodies"], c):
            body.set_facecolor(cc)
            body.set_alpha(0.45)
        ax.boxplot(data, positions=range(len(sigmas)), widths=0.25, showfliers=False,
                   medianprops={"color": "black"})
        for i, d in enumerate(data):
            ax.scatter(np.full(len(d), i) + np.linspace(-0.12, 0.12, len(d)), d, color=c[i], edgecolor="black",
                       s=28, zorder=3)
        ax.set_xticks(range(len(sigmas)))
        ax.set_xticklabels([f"{s:g}" for s in sigmas])
        ax.set_xlabel("log-normal spread $\\sigma$ of the noise-model error")
        ax.set_ylabel(lab)
    axes[0].set_title("Decision errors (5 seeds, points)", fontsize=10)
    axes[1].set_title("Cost of those errors", fontsize=10)
    save(fig_, fig / "fig5_calibration_regret.png")
    g = df.groupby("sigma")[["misclassified", "mean_regret"]].agg(["mean", "std"])
    return {f"sigma={s:g}": {"misclassified_mean": float(g.loc[s, ("misclassified", "mean")]),
                             "misclassified_std": float(g.loc[s, ("misclassified", "std")]),
                             "regret_mean": float(g.loc[s, ("mean_regret", "mean")]),
                             "regret_std": float(g.loc[s, ("mean_regret", "std")])} for s in sigmas}


def exp_headline(res: Path) -> dict:
    """Headline table at the default operating point and vs the unencoded qubit."""
    rows = []
    for tau in (2.0, 10.0, 30.0, 100.0):
        prm = BASE.with_(tau_m=tau)
        rows.append({"tau_m": tau, "mb": measurement_based_fail(prm), "mf": measurement_free_fail(prm),
                     "bare_mb_duration": bare_fail(prm, total_time_mb(prm)),
                     "bare_mf_duration": bare_fail(prm, total_time_mf(prm)),
                     "tau_star": tau_star(prm)})
    pd.DataFrame(rows).to_csv(res / "headline.csv", index=False)
    # pseudo-threshold vs unencoded qubit
    pt = {}
    for scheme, fn, tf in (("mb", measurement_based_fail, total_time_mb), ("mf", measurement_free_fail, total_time_mf)):
        lo, hi = 1e-6, 5e-2
        prm = BASE.with_(tau_m=10.0)
        for _ in range(40):
            mid = np.sqrt(lo * hi)
            q = prm.with_(p2=mid)
            if fn(q) < bare_fail(q, tf(q)):
                lo = mid
            else:
                hi = mid
        pt[scheme] = float(np.sqrt(lo * hi))
    return {"pseudo_threshold_p2_vs_unencoded_tau10": pt}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="< 1 minute smoke run")
    args = ap.parse_args()
    res, fig = ROOT / "results", ROOT / "figures"
    if args.quick:  # smoke runs never overwrite the published full results
        res, fig = res / "quick", fig / "quick"
    res.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    summary: dict = {"quick": args.quick, "base_params": BASE.__dict__ | {}}
    summary["validation"] = exp_validation(res, args.quick); timed("validation", t0)
    summary["headline"] = exp_headline(res); timed("headline", t0)
    summary["phase_maps"] = exp_phase_maps(res, fig, args.quick); timed("phase maps", t0)
    summary["tau_star"] = exp_tau_star(res, fig, args.quick); timed("break-even sweeps", t0)
    summary["finite_shots"] = exp_finite_shots(res, fig, args.quick); timed("finite shots", t0)
    summary["rounds"] = exp_rounds(res, fig, args.quick); timed("cycle scaling", t0)
    summary["readout_threshold"] = exp_thresholds(res); timed("readout thresholds", t0)
    summary["calibration"] = exp_calibration(res, fig, args.quick); timed("calibration regret", t0)
    summary["runtime_seconds"] = time.time() - t0
    out = res / "summary.json"
    out.write_text(json.dumps(summary, indent=2, default=float))
    print(f"done in {time.time() - t0:.1f}s -> {out}")


if __name__ == "__main__":
    main()
