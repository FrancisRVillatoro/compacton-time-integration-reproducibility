"""Figure 3: controlled damping with the generalized-alpha family (p = 3,
c1 = 3, c2 = -1, N = 1000, dt = 0.2, mu = 1e-7, window [9, 21]).

Left  (fig3_controlled_damping_tradeoff_left.pdf):  high-wavenumber AUC
      versus relative drift of I2; generalized-alpha points joined in
      decreasing order of rho_inf and labelled by its value; Pade method.
Right (fig3_controlled_damping_tradeoff_right.pdf): the same quantities
      divided by their Pade values, with the unit lines.

Data: data/fig3_controlled_damping_fine_corrected.csv.
"""
from __future__ import annotations

import numpy as np

from figure_style import GA_COLOR, GREY, METHOD, log_axis, log_axis_short, panel, rows, save

ORDER = (0.9, 0.75, 0.5, 0.25, 0.1, 0.0)


def load():
    raw = rows("fig3_controlled_damping_fine_corrected.csv")
    ga = {float(r["rho_inf"]): r for r in raw
          if r["label"].startswith("GA_") and r["status"] == "ok"}
    pade = next(r for r in raw if r["label"] == "M2")
    drift = np.array([float(ga[r]["abs_I2_drift"]) for r in ORDER])
    auc = np.array([float(ga[r]["auc"]) for r in ORDER])
    return drift, auc, float(pade["abs_I2_drift"]), float(pade["auc"])


def draw(ax, x, y, x_pade, y_pade, offsets, pade_offset):
    ax.plot(x, y, "-", color=GA_COLOR, lw=0.9, zorder=2)
    ax.plot(x, y, "o", color=GA_COLOR, ms=4.2, mfc="white", mew=0.9, zorder=3,
            label=r"generalized-$\alpha$ ($\rho_\infty$ labelled)")
    for rho, xi, yi in zip(ORDER, x, y):
        dx, dy, ha, va = offsets[rho]
        ax.annotate(rf"${rho:g}$", (xi, yi), xytext=(dx, dy),
                    textcoords="offset points", ha=ha, va=va, fontsize=8.0,
                    color=GA_COLOR)
    st = METHOD["pade"]
    ax.plot([x_pade], [y_pade], ls="none", marker=st["marker"], color=st["color"],
            ms=5.0, zorder=5, label=st["label"])
    dx, dy, ha, va = pade_offset
    ax.annotate(st["label"], (x_pade, y_pade), xytext=(dx, dy),
                textcoords="offset points", ha=ha, va=va, fontsize=8.0,
                color=st["color"])


def main():
    drift, auc, pade_drift, pade_auc = load()

    # ------------------------------------------------------------ left panel
    fig, ax = panel()
    offsets = {0.9: (-5, -2, "right", "top"), 0.75: (-5, 2, "right", "bottom"),
               0.5: (-5, 2, "right", "bottom"), 0.25: (3, 5, "left", "bottom"),
               0.1: (3, 5, "left", "bottom"), 0.0: (3, 5, "left", "bottom")}
    draw(ax, drift, auc, pade_drift, pade_auc, offsets, (0, -7, "center", "top"))
    ax.set_xscale("log")
    ax.set_yscale("log")
    log_axis(ax, "x")
    log_axis(ax, "y")
    ax.set_xlim(5.0e-4, 8.0e-2)
    ax.set_ylim(1.0e-4, 1.0e-1)
    ax.set_xlabel(r"$|\Delta I_2/I_2|$")
    ax.set_ylabel(r"$\rho_{\mathrm{hw}}^{\mathrm{AUC}}$")
    ax.legend(loc="upper right")
    ax.text(0.03, 0.04, "lower is better", transform=ax.transAxes, fontsize=8.0,
            color=GREY)
    save(fig, "fig3_controlled_damping_tradeoff_left")

    # ----------------------------------------------------------- right panel
    fig, ax = panel()
    x = drift / pade_drift
    y = auc / pade_auc
    ax.axhline(1.0, ls="--", color=GREY, lw=0.8, zorder=1)
    ax.axvline(1.0, ls="--", color=GREY, lw=0.8, zorder=1)
    offsets = {0.9: (-5, -2, "right", "top"), 0.75: (-5, 2, "right", "bottom"),
               0.5: (-5, 2, "right", "bottom"), 0.25: (3, 5, "left", "bottom"),
               0.1: (3, 5, "left", "bottom"), 0.0: (3, 5, "left", "bottom")}
    draw(ax, x, y, 1.0, 1.0, offsets, (-5, -7, "right", "top"))
    ax.set_xscale("log")
    ax.set_yscale("log")
    log_axis_short(ax, "x", subs=(1.0, 2.0, 5.0))
    log_axis(ax, "y")
    ax.set_xlim(0.15, 20.0)
    ax.set_ylim(0.15, 150.0)
    ax.set_xlabel("drift ratio to Padé")
    ax.set_ylabel("AUC ratio to Padé")
    ax.legend(loc="upper right")
    ax.text(0.03, 0.04, "lower is better", transform=ax.transAxes, fontsize=8.0,
            color=GREY)
    save(fig, "fig3_controlled_damping_tradeoff_right")


if __name__ == "__main__":
    main()
