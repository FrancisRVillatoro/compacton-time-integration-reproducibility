"""Figure 2: mesh scaling of the stiffest Jacobian mode (p = 3).

Left  (fig2_jacobian_h3_scaling_left.pdf):  |lambda_max| versus h for the
      separated and the overlapping compacton-anticompacton pair, with
      power laws fitted on the three finest grids and an h^-3 guide.
Right (fig2_jacobian_h3_scaling_right.pdf): amplification of the stiffest
      mode of the separated pair for five integrators, with h^3 and h^(3/2)
      guides.

Data: data/fig2_jacobian_h3_scaling.csv.
"""
from __future__ import annotations

import numpy as np

from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

from figure_style import GREY, METHOD, column, guide, log_axis, panel, rows, save


def h_ticks(ax, h):
    """Tick marks at the five mesh sizes only."""
    ax.xaxis.set_major_locator(FixedLocator(list(h)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())


def main():
    data = rows("fig2_jacobian_h3_scaling.csv")
    sep = sorted((r for r in data if r["state"] == "separated"), key=lambda r: float(r["h"]))
    ov = sorted((r for r in data if r["state"] == "overlap"), key=lambda r: float(r["h"]))
    h = column(sep, "h")
    lam_s = column(sep, "abs_lambda_max")
    lam_o = column(ov, "abs_lambda_max")
    # the three finest grids define the fitted exponents quoted in the text
    slope_s = -np.polyfit(np.log(h[:3]), np.log(lam_s[:3]), 1)[0]
    slope_o = -np.polyfit(np.log(h[:3]), np.log(lam_o[:3]), 1)[0]

    # ------------------------------------------------------------ left panel
    fig, ax = panel()
    ax.loglog(h, lam_s, "-", marker="s", color=METHOD["pade"]["color"], mfc="white",
              label=rf"separated pair, fit $h^{{-{slope_s:.2f}}}$")
    ax.loglog(h, lam_o, "-", marker="D", color=METHOD["midpoint"]["color"], mfc="white",
              label=rf"overlapping pair, fit $h^{{-{slope_o:.2f}}}$")
    guide(ax, h, 0.45 * lam_s[0], -3.0, r"$\propto h^{-3}$", x_anchor=h[0],
          x_text=0.10, y_text_factor=0.30, ha="center", va="top")
    h_ticks(ax, h)
    log_axis(ax, "y")
    ax.set_xlim(0.021, 0.48)
    ax.set_ylim(3.0e2, 2.0e7)
    ax.set_xlabel(r"mesh size $h$")
    ax.set_ylabel(r"$|\lambda_{\max}|$")
    ax.legend(loc="upper right", handlelength=1.8)
    save(fig, "fig2_jacobian_h3_scaling_left")

    # ----------------------------------------------------------- right panel
    # Curves are labelled directly, next to the curve, instead of a legend.
    fig, ax = panel()
    series = [("M1_amp", "midpoint", None, None),
              ("M3_amp", "composition", None, None),
              ("GA_rho05_amp", "ga05", 0.10, (0.55, "top")),
              ("GA_rho0_amp", "ga0", 0.10, (1.9, "bottom")),
              ("M2_amp", "pade", 0.20, (0.18, "top"))]
    for key, name, xt, place in series:
        st = METHOD[name]
        y = column(sep, key)
        ax.loglog(h, y, "-", marker=st["marker"], color=st["color"],
                  mfc="white" if name == "composition" else st["color"], zorder=3)
        if xt is not None:
            i = int(np.argmin(np.abs(h - xt)))
            ax.text(h[i], y[i] * place[0], st["label"], color=st["color"], fontsize=8.2,
                    ha="center", va=place[1])
    i = int(np.argmin(np.abs(h - 0.10)))
    ax.text(h[i], 1.7, "implicit midpoint, composition", color=METHOD["midpoint"]["color"],
            fontsize=8.2, ha="center", va="bottom")
    guide(ax, h, 0.45 * column(sep, "M2_amp")[0], 3.0, r"$\propto h^{3}$",
          x_anchor=h[0], x_text=0.30, y_text_factor=0.25, ha="center", va="top")
    guide(ax, h, 0.45 * column(sep, "GA_rho0_amp")[0], 1.5, r"$\propto h^{3/2}$",
          x_anchor=h[0], x_text=0.042, y_text_factor=0.40, ha="center", va="top")
    h_ticks(ax, h)
    log_axis(ax, "y")
    ax.set_xlim(0.021, 0.48)
    ax.set_ylim(5.0e-6, 8.0)
    ax.set_xlabel(r"mesh size $h$")
    ax.set_ylabel("stiff-mode amplification")
    save(fig, "fig2_jacobian_h3_scaling_right")


if __name__ == "__main__":
    main()
