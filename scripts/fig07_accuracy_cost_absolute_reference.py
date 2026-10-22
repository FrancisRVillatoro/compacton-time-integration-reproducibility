"""Accuracy-cost figure (file fig7_accuracy_cost_absolute_reference.pdf):
coherent-core error against the N = 3200 absolute-flux Direct-WENO
reference versus wall time for the nine N = 1600 computations of the
L-stable comparison, joined per method in decreasing order of the time
step, with the resolution indicator of the reference as a dashed line.

Data: data/table3_final.csv.
"""
from __future__ import annotations

from figure_style import GREY, METHOD, log_axis, panel, rows, save

RESOLUTION_INDICATOR = 4.426e-4


def main():
    data = rows("table3_final.csv")
    fig, ax = panel(width_mm=110.0, height_mm=70.0, left=15.0)
    for key, name in (("M2", "pade"), ("RadauIIA2", "radau"), ("TR-BDF2", "trbdf2")):
        st = METHOD[name]
        rr = sorted((r for r in data if r["method"] == key),
                    key=lambda r: float(r["dt"]), reverse=True)
        x = [float(r["wall_s"]) for r in rr]
        y = [float(r["E_core"]) for r in rr]
        ax.plot(x, y, "-", marker=st["marker"], color=st["color"], label=st["label"],
                zorder=3)
        # only the largest step of each method is labelled; the table gives the rest
        dt0 = float(rr[0]["dt"])
        ax.annotate(rf"$\Delta t={dt0:g}$", (x[0], y[0]),
                    xytext=(6, 0) if name == "radau" else (0, 6),
                    textcoords="offset points",
                    ha="left" if name == "radau" else "center",
                    va="center" if name == "radau" else "bottom",
                    fontsize=7.8, color=st["color"])
    ax.axhline(RESOLUTION_INDICATOR, ls="--", color=GREY, lw=0.8, zorder=1)
    ax.text(15.6, RESOLUTION_INDICATOR * 1.15, "resolution indicator of the reference",
            ha="right", va="bottom", fontsize=8.0, color=GREY)
    ax.set_yscale("log")
    log_axis(ax, "y")
    ax.set_xlim(2.0, 16.0)
    ax.set_ylim(3.0e-4, 4.0e-2)
    ax.set_xlabel("wall time [s]")
    ax.set_ylabel(r"$E_{\mathrm{core}}$")
    ax.legend(loc="upper left", handlelength=1.8)
    save(fig, "fig7_accuracy_cost_absolute_reference")


if __name__ == "__main__":
    main()
