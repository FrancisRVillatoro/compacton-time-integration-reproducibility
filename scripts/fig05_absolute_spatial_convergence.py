"""Figure 5: Direct-WENO spatial validation (p = 2, absolute flux).

Left  (fig5_absolute_spatial_convergence_left.pdf):  relative L2 errors of
      isolated compactons with c = 3 and c = 1 versus h, with an h^2 guide.
Right (fig5_absolute_spatial_convergence_right.pdf): pairwise differences
      between successive grids of the same-sign collision at T = 10
      (registered field, modes |m| <= 32, fitted cores) versus the number
      of nodes of the finer grid.

Data: data/isolated_compacton_validation.csv,
      data/absolute_spatial_convergence.csv.
"""
from __future__ import annotations

import numpy as np

from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

from figure_style import BLUE, GREEN, VERMILION, column, guide, log_axis, panel, rows, save


def main():
    iso = rows("isolated_compacton_validation.csv")
    con = rows("absolute_spatial_convergence.csv")

    # ------------------------------------------------------------ left panel
    h = column(iso, "h")
    e3 = column(iso, "relL2_c3")
    e1 = column(iso, "relL2_c1")
    fig, ax = panel()
    ax.loglog(h, e3, "-", marker="o", color=BLUE, label=r"$c=3$")
    ax.loglog(h, e1, "-", marker="s", color=VERMILION, label=r"$c=1$")
    guide(ax, h, 0.55 * e1[-1], 2.0, r"$\propto h^{2}$", x_anchor=h[-1],
          x_text=0.11, y_text_factor=0.55, ha="center")
    ax.xaxis.set_major_locator(FixedLocator(list(h)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    log_axis(ax, "y")
    ax.set_xlim(0.05, 0.32)
    ax.set_ylim(2.0e-5, 3.0e-3)
    ax.set_xlabel(r"mesh size $h$")
    ax.set_ylabel(r"relative $L^2$ error")
    ax.legend(loc="upper left", handlelength=1.8)
    save(fig, "fig5_absolute_spatial_convergence_left")

    # ----------------------------------------------------------- right panel
    nf = column(con, "N_fine")
    ef = column(con, "registered_relL2")
    el = column(con, "low32_registered")
    ec = column(con, "E_core")
    p_field = np.log2(ef[0] / ef[1])
    p_low = np.log2(el[0] / el[1])
    p_core = np.log2(ec[0] / ec[1])
    fig, ax = panel()
    ax.semilogy(nf, ef, "-", marker="o", color=BLUE,
                label=rf"registered field, order ${p_field:.2f}$")
    ax.semilogy(nf, el, "-", marker="s", color=VERMILION,
                label=rf"modes $|m|\leq32$, order ${p_low:.2f}$")
    ax.semilogy(nf, ec, "-", marker="^", color=GREEN,
                label=rf"fitted cores, order ${p_core:.2f}$")
    log_axis(ax, "y")
    ax.set_xticks([1600, 3200])
    ax.set_xticklabels(["1600", "3200"])
    ax.tick_params(axis="x", which="minor", bottom=False, top=False)
    ax.set_xlim(1300, 3500)
    ax.set_ylim(2.0e-4, 2.0e-2)
    ax.set_xlabel(r"nodes of the finer grid, $N$")
    ax.set_ylabel("pairwise relative difference")
    ax.legend(loc="upper right", handlelength=1.8)
    save(fig, "fig5_absolute_spatial_convergence_right")


if __name__ == "__main__":
    main()
