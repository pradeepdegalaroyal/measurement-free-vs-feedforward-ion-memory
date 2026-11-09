"""Figure helpers (viridis palette, 200 dpi)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

DPI = 200
CMAP = plt.get_cmap("viridis")


def colors(n: int) -> list[tuple[float, float, float, float]]:
    """``n`` distinct viridis colours (dark end to light end, skipping the pale yellow)."""
    return [CMAP(v) for v in np.linspace(0.05, 0.85, n)]


def save(fig: plt.Figure, path: Path) -> None:
    """Tight-layout and write a 200 dpi PNG."""
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def phase_panels(maps: dict[int, np.ndarray], p2: np.ndarray, tau: np.ndarray, path: Path) -> None:
    """Filled-contour phase diagrams of ``log10(P_MB / P_MF)``; zero contour is the boundary."""
    n = len(maps)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 4.2), sharey=True)
    axes = np.atleast_1d(axes)
    allv = np.concatenate([m.ravel() for m in maps.values()])
    lim = float(np.ceil(max(abs(allv.min()), abs(allv.max())) * 4) / 4)
    levels = np.linspace(-lim, lim, 17)
    cf = None
    for ax, (n2, m) in zip(axes, maps.items()):
        cf = ax.contourf(p2, tau, m, levels=levels, cmap="viridis", extend="both")
        ax.contour(p2, tau, m, levels=[0.0], colors="white", linewidths=2.2, linestyles="--")
        ax.contour(p2, tau, m, levels=[0.0], colors="black", linewidths=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{n2} two-qubit gate{'s' if n2 > 1 else ''} per Toffoli")
        ax.set_xlabel("two-qubit error probability $p_2$")
    axes[0].set_ylabel("measurement time $\\tau_m$ (two-qubit gate-times)")
    cb = fig.colorbar(cf, ax=list(axes), shrink=0.9, pad=0.02)
    cb.set_label("$\\log_{10}(P_{\\mathrm{feed\\text{-}forward}} / P_{\\mathrm{measurement\\text{-}free}})$"
                 "  (> 0: measurement-free better)")
    fig.suptitle("Where coherent correction beats feed-forward (3 cycles, exact)", y=1.02)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
