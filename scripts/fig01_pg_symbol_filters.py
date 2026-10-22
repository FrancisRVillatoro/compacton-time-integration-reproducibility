"""Figure 1: exact Petrov-Galerkin symbol and amplification profiles.

Left  (fig1_pg_symbol_filters_left.pdf):  Im z(theta) for the exact symbol
      with q = 5.625 and q = 0, and for the continuum proxy with q = 5.625.
Right (fig1_pg_symbol_filters_right.pdf): |R(z(theta))| along the exact
      locus with q = 5.625 for the six integrators, and the Pade method on
      the continuum proxy.

Data: data/fig1_pg_symbol_filters.csv (h = 0.05, dt = 0.1, mu = 1e-7, c0 = 1).
The generalized-alpha amplification is the spectral radius of the 2 x 2
amplification matrix (generalized_alpha.py).
"""
from __future__ import annotations

import numpy as np

from figure_style import GREY, METHOD, column, log_axis, panel, rows, save
from generalized_alpha import ga_rho


def main():
    data = rows("fig1_pg_symbol_filters.csv")
    x = column(data, "theta_over_pi")
    z_exact = column(data, "z_q5625_real") + 1j * column(data, "z_q5625_imag")
    im_exact = z_exact.imag / 1.0e4
    im_cont = column(data, "z_cont_imag") / 1.0e4
    im_q0 = column(data, "z_q0_imag") / 1.0e4

    # ------------------------------------------------------------ left panel
    fig, ax = panel()
    ax.plot(x, im_exact, "-", color=METHOD["midpoint"]["color"],
            label=r"exact symbol, $q=5.625$")
    ax.plot(x, im_cont, "--", color=METHOD["pade"]["color"],
            label=r"continuum proxy, $q=5.625$")
    ax.plot(x, im_q0, ":", color=GREY, label=r"exact symbol, $q=0$")
    i = int(np.argmax(im_exact))
    ax.plot(x[i], im_exact[i], "o", ms=4.5, mfc="white",
            mec=METHOD["midpoint"]["color"], zorder=5)
    ax.annotate(r"$\theta^*\simeq0.80\pi$", xy=(x[i], im_exact[i]),
                xytext=(0.60, 7.2), fontsize=8.5,
                arrowprops=dict(arrowstyle="-", lw=0.6, color=GREY))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.5, 14.5)
    ax.set_xlabel(r"$\theta/\pi$")
    ax.set_ylabel(r"$\mathrm{Im}\,z(\theta)/10^{4}$")
    ax.legend(loc="upper left")
    save(fig, "fig1_pg_symbol_filters_left")

    # ----------------------------------------------------------- right panel
    # Curves are labelled directly, next to the curve, instead of a legend.
    fig, ax = panel()
    ga05 = np.array([ga_rho(complex(v), 0.5) for v in z_exact])
    ga0 = np.array([ga_rho(complex(v), 0.0) for v in z_exact])
    curves = [
        (column(data, "amp_midpoint_exact"), "-", METHOD["midpoint"]["color"],
         "implicit midpoint, composition", 0.50, 1.6, "center", "bottom"),
        (ga05, "-", METHOD["ga05"]["color"], METHOD["ga05"]["label"],
         0.50, 0.55, "center", "top"),
        (ga0, "-", METHOD["ga0"]["color"], METHOD["ga0"]["label"],
         0.50, 1.9, "center", "bottom"),
        (column(data, "amp_TRBDF2_exact"), "-", METHOD["trbdf2"]["color"], "TR-BDF2",
         0.42, 1.9, "center", "bottom"),
        (column(data, "amp_M2_Radau_exact"), "-", METHOD["pade"]["color"],
         METHOD["pade"]["label"] + ", Radau IIA", 0.38, 0.30, "right", "top"),
        (column(data, "amp_M2_continuum"), "--", GREY,
         METHOD["pade"]["label"] + ", continuum proxy", 0.72, 0.45, "center", "top"),
    ]
    for y, ls, color, label, xt, factor, ha, va in curves:
        ax.semilogy(x, y, ls, color=color)
        i = int(np.argmin(np.abs(x - xt)))
        ax.text(x[i], y[i] * factor, label, color=color, fontsize=8.2, ha=ha, va=va)
    log_axis(ax, "y")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(3.0e-6, 8.0)
    ax.set_xlabel(r"$\theta/\pi$")
    ax.set_ylabel(r"amplification $|R(z(\theta))|$")
    save(fig, "fig1_pg_symbol_filters_right")


if __name__ == "__main__":
    main()
