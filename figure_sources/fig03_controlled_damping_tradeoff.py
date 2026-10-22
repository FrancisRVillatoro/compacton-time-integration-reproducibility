"""Figure 3: corrected fine-grid generalized-alpha damping/accuracy tradeoff."""
from __future__ import annotations

import numpy as np

from plot_support import output, rows, single_panel


ORDER = (0.9, 0.75, 0.5, 0.25, 0.1, 0.0)


def load_data():
    raw = rows("fig3_controlled_damping_fine_corrected.csv")
    ga_by_rho = {
        float(row["rho_inf"]): row
        for row in raw
        if row["label"].startswith("GA_") and row["status"] == "ok"
    }
    ga = [ga_by_rho[rho] for rho in ORDER]
    pade = next(row for row in raw if row["label"] == "M2")
    return ga, pade


def annotate_rho(ax, x, y, normalized=False):
    if normalized:
        offsets = {
            0.9: (5, -16),
            0.75: (-18, 6),
            0.5: (-18, 6),
            0.25: (7, -12),
            0.1: (7, 5),
            0.0: (7, 5),
        }
    else:
        offsets = {
            0.9: (5, -16),
            0.75: (-18, 6),
            0.5: (-18, 7),
            0.25: (7, 7),
            0.1: (7, 7),
            0.0: (7, 4),
        }
    for rho, xx, yy in zip(ORDER, x, y):
        ax.annotate(
            fr"${rho:g}$",
            (xx, yy),
            xytext=offsets[rho],
            textcoords="offset points",
            fontsize=7.9,
        )


def main():
    ga, pade = load_data()
    drift = np.array([float(row["abs_I2_drift"]) for row in ga])
    auc = np.array([float(row["auc"]) for row in ga])
    pade_drift = float(pade["abs_I2_drift"])
    pade_auc = float(pade["auc"])

    # Both axes are error diagnostics; the lower-left direction is favorable.
    fig, ax = single_panel()
    ax.loglog(
        drift,
        auc,
        "o-",
        mfc="white",
        color="C0",
        label=r"generalized-$\alpha$ ($\rho_\infty$ labels)",
        markersize=5.2,
    )
    ax.loglog(
        [pade_drift],
        [pade_auc],
        "s",
        color="C1",
        label=r"$R_{1,2}$ Pad\'e",
        markersize=6.8,
    )
    annotate_rho(ax, drift, auc, normalized=False)
    ax.annotate(
        r"$R_{1,2}$ Pad\'e",
        (pade_drift, pade_auc),
        xytext=(7, -16),
        textcoords="offset points",
        fontsize=8.0,
    )
    ax.set_xlim(7.0e-4, 6.0e-2)
    ax.set_ylim(1.4e-4, 7.0e-2)
    ax.set_xlabel(r"$|\Delta I_2/I_2|$")
    ax.set_ylabel(r"post-collision $\rho_{\rm hw}^{\rm post,AUC}$")
    ax.legend(frameon=False, loc="upper right", fontsize=8.2, handlelength=1.4)
    ax.text(
        0.04,
        0.045,
        "lower is better",
        transform=ax.transAxes,
        fontsize=7.5,
        color="0.35",
    )
    output(fig, "fig3_controlled_damping_tradeoff_left.pdf")

    fig, ax = single_panel()
    x = drift / pade_drift
    y = auc / pade_auc
    ax.axhline(1.0, ls="--", color="0.62", lw=1.0)
    ax.axvline(1.0, ls="--", color="0.62", lw=1.0)
    ax.loglog(
        x,
        y,
        "o-",
        mfc="white",
        color="C0",
        label=r"generalized-$\alpha$ ($\rho_\infty$ labels)",
        markersize=5.2,
    )
    ax.loglog(
        [1.0],
        [1.0],
        "s",
        color="C1",
        label=r"$R_{1,2}$ Pad\'e",
        markersize=6.8,
    )
    annotate_rho(ax, x, y, normalized=True)
    ax.annotate(
        r"$R_{1,2}$ Pad\'e",
        (1.0, 1.0),
        xytext=(7, 7),
        textcoords="offset points",
        fontsize=8.0,
    )
    ax.set_xlim(0.18, 15.0)
    ax.set_ylim(0.18, 105.0)
    ax.set_xlabel(r"invariant-drift ratio to $R_{1,2}$ Pad\'e")
    ax.set_ylabel(r"AUC ratio to $R_{1,2}$ Pad\'e")
    ax.legend(frameon=False, loc="upper right", fontsize=8.2, handlelength=1.4)
    ax.text(
        0.04,
        0.045,
        "lower is better",
        transform=ax.transAxes,
        fontsize=7.5,
        color="0.35",
    )
    output(fig, "fig3_controlled_damping_tradeoff_right.pdf")


if __name__ == "__main__":
    main()
