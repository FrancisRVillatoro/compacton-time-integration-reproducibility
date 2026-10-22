"""Figure 4: high-throughput robustness summary.

Left  (fig4_highthroughput_summary_left.pdf):  fraction of completed runs
      per method in the primary, collision-severity and long-time banks.
Right (fig4_highthroughput_summary_right.pdf): median (filled) and 95th
      percentile (open) of the high-wavenumber AUC over completed runs of
      the primary and collision-severity banks.

Data: data/fig4_highthroughput_summary.csv.
"""
from __future__ import annotations

from matplotlib.lines import Line2D
import numpy as np

from figure_style import METHOD, log_axis, panel, rows, save

METHODS = (("M1", "midpoint"), ("M2", "pade"), ("M3", "composition"))
BANKS = ("mu-dt", "severity", "long-time")
EN_DASH = "–"
BANK_LABELS = {"mu-dt": r"primary $\mu$" + EN_DASH + r"$\Delta t$",
               "severity": "collision severity",
               "long-time": "long time"}


def load():
    return {(r["bank"], r["method"]): r for r in rows("fig4_highthroughput_summary.csv")}


def main():
    data = load()

    # ------------------------------------------------------------ left panel
    fig, ax = panel()
    centers = np.arange(len(BANKS), dtype=float)
    width = 0.25
    for j, (key, name) in enumerate(METHODS):
        st = METHOD[name]
        frac = np.array([float(data[(b, key)]["completed"]) / float(data[(b, key)]["assigned"])
                         for b in BANKS])
        xpos = centers + (j - 1) * width
        bars = ax.bar(xpos, frac, width=width * 0.92, color=st["color"], label=st["label"],
                      zorder=3)
        for bar, value in zip(bars, frac):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 0.02, f"{100 * value:.0f}%",
                    ha="center", va="bottom", fontsize=7.8)
    ax.set_xticks(centers)
    ax.set_xticklabels([BANK_LABELS[b] for b in BANKS])
    ax.tick_params(axis="x", which="both", top=False)
    ax.tick_params(axis="x", which="minor", bottom=False)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(0.0, 1.36)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylabel("fraction of completed runs")
    ax.legend(loc="upper center", ncol=3, columnspacing=0.9, handlelength=1.2,
              handletextpad=0.4, bbox_to_anchor=(0.5, 1.0))
    save(fig, "fig4_highthroughput_summary_left")

    # ----------------------------------------------------------- right panel
    fig, ax = panel()
    auc_banks = BANKS[:2]
    base = np.arange(len(auc_banks), dtype=float)
    offsets = {"M1": -0.22, "M2": 0.0, "M3": 0.22}
    for key, name in METHODS:
        st = METHOD[name]
        med = np.array([float(data[(b, key)]["median_auc"]) for b in auc_banks])
        p95 = np.array([float(data[(b, key)]["p95_auc"]) for b in auc_banks])
        xpos = base + offsets[key]
        for xx, m, u in zip(xpos, med, p95):
            ax.plot([xx, xx], [m, u], "-", color=st["color"], lw=1.0, zorder=2)
        ax.plot(xpos, med, ls="none", marker=st["marker"], color=st["color"], ms=5.0, zorder=3)
        ax.plot(xpos, p95, ls="none", marker=st["marker"], mfc="white", mec=st["color"],
                color=st["color"], ms=5.0, zorder=3)
    ax.set_yscale("log")
    log_axis(ax, "y")
    ax.set_xticks(base)
    ax.set_xticklabels([BANK_LABELS[b] for b in auc_banks])
    ax.tick_params(axis="x", which="both", top=False)
    ax.tick_params(axis="x", which="minor", bottom=False)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(1.0e-6, 1.0e3)
    ax.set_ylabel(r"$\rho_{\mathrm{hw}}^{\mathrm{AUC}}$")
    method_handles = [Line2D([0], [0], ls="none", marker=METHOD[n]["marker"],
                             color=METHOD[n]["color"], ms=5.0, label=METHOD[n]["label"])
                      for _, n in METHODS]
    value_handles = [Line2D([0], [0], ls="none", marker="o", color="black", ms=5.0,
                            label="median"),
                     Line2D([0], [0], ls="none", marker="o", mfc="white", mec="black",
                            ms=5.0, label="95th percentile")]
    first = ax.legend(handles=method_handles, loc="upper left", handlelength=1.2)
    ax.add_artist(first)
    ax.legend(handles=value_handles, loc="upper right", handlelength=1.2)
    save(fig, "fig4_highthroughput_summary_right")


if __name__ == "__main__":
    main()
