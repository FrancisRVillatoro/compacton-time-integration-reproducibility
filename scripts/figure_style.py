"""Shared style for the manuscript figures.

Every panel is a separate PDF of identical size (88 mm x 64 mm, one journal
column) with the same fonts, line widths and method-to-colour assignment, so
that any two panels placed side by side at 0.49\\linewidth align.  Only fonts
bundled with matplotlib are used (STIX), so the scripts reproduce identically
on any machine without a TeX installation.  No panel carries a printed title
or an (a)/(b) label; panels are identified by position in the captions.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (FuncFormatter, LogLocator, LogFormatterSciNotation,
                               NullFormatter)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"

# ---------------------------------------------------------------- geometry
MM = 1.0 / 25.4
PANEL_W = 88.0 * MM
PANEL_H = 64.0 * MM
MARGIN = dict(left=14.5, bottom=11.0, right=2.0, top=2.0)   # mm

# ---------------------------------------------------------------- palette
# Okabe-Ito colour-blind-safe palette
BLACK = "#000000"
ORANGE = "#E69F00"
SKY = "#56B4E9"
GREEN = "#009E73"
BLUE = "#0072B2"
VERMILION = "#D55E00"
PURPLE = "#CC79A7"
GREY = "#767676"

# one visual identity per method, used in every figure; the labels are the
# names used in the text of the manuscript
METHOD = {
    "midpoint":    dict(color=BLUE,      marker="o", label="implicit midpoint"),
    "pade":        dict(color=VERMILION, marker="s", label="Padé $R_{1,2}$"),
    "composition": dict(color=SKY,       marker="^", label="composition"),
    "radau":       dict(color=GREEN,     marker="D", label="Radau IIA"),
    "trbdf2":      dict(color=PURPLE,    marker="v", label="TR-BDF2"),
    "ga0":         dict(color=ORANGE,    marker="p", label=r"generalized-$\alpha$, $\rho_\infty=0$"),
    "ga05":        dict(color=BLACK,     marker="h", label=r"generalized-$\alpha$, $\rho_\infty=0.5$"),
}
GA_COLOR = ORANGE

mpl.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 9.0,
    "axes.labelsize": 10.0,
    "axes.titlesize": 10.0,
    "xtick.labelsize": 9.0,
    "ytick.labelsize": 9.0,
    "legend.fontsize": 8.5,
    "legend.frameon": False,
    "legend.handlelength": 2.2,
    "legend.handletextpad": 0.6,
    "legend.labelspacing": 0.35,
    "legend.borderaxespad": 0.4,
    "axes.linewidth": 0.6,
    "axes.formatter.use_mathtext": True,
    "axes.formatter.limits": (-3, 4),
    "lines.linewidth": 1.1,
    "lines.markersize": 4.2,
    "lines.markeredgewidth": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.2,
    "ytick.major.size": 3.2,
    "xtick.minor.size": 1.8,
    "ytick.minor.size": 1.8,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.45,
    "ytick.minor.width": 0.45,
    "xtick.major.pad": 3.0,
    "ytick.major.pad": 3.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "standard",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


def rows(name: str) -> list[dict[str, str]]:
    """Read data/<name> as a list of dictionaries (one per row)."""
    with (DATA / name).open("r", newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def column(records: list[dict[str, str]], field: str) -> np.ndarray:
    return np.array([float(r[field]) for r in records], dtype=float)


def panel(width_mm: float = 88.0, height_mm: float = 64.0, **margins):
    """Return (fig, ax) for one panel with a fixed axes box."""
    m = dict(MARGIN)
    m.update(margins)
    w, h = width_mm * MM, height_mm * MM
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([m["left"] * MM / w,
                       m["bottom"] * MM / h,
                       1.0 - (m["left"] + m["right"]) * MM / w,
                       1.0 - (m["bottom"] + m["top"]) * MM / h])
    ax.minorticks_on()
    return fig, ax


def log_axis(ax, axis="y"):
    """Clean powers-of-ten ticks on a logarithmic axis."""
    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_major_locator(LogLocator(base=10.0, numticks=12))
    a.set_major_formatter(LogFormatterSciNotation(base=10.0))
    a.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=12))
    a.set_minor_formatter(NullFormatter())


def log_axis_short(ax, axis="x", subs=(1.0, 2.0, 5.0)):
    """Decimal labels at 1, 2, 5 x 10^k for logarithmic axes spanning
    fewer than two decades."""
    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_major_locator(LogLocator(base=10.0, subs=subs, numticks=20))
    a.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    a.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=20))
    a.set_minor_formatter(NullFormatter())


def guide(ax, x, y_anchor, exponent, label, x_anchor=None, x_text=None,
          y_text_factor=1.0, color=GREY, ha="center", va="top"):
    """Dashed reference line y = C x^exponent through (x_anchor, y_anchor),
    labelled at x_text with a multiplicative vertical offset."""
    x = np.asarray(x, dtype=float)
    if x_anchor is None:
        x_anchor = x[-1]
    c = y_anchor / x_anchor ** exponent
    xx = np.array([x.min(), x.max()])
    ax.plot(xx, c * xx ** exponent, "--", color=color, lw=0.9, zorder=1)
    if x_text is None:
        x_text = np.sqrt(xx[0] * xx[1])
    ax.text(x_text, c * x_text ** exponent * y_text_factor, label, color=color,
            fontsize=8.5, ha=ha, va=va)


def save(fig, stem: str) -> Path:
    """Write figures/<stem>.pdf with fixed metadata and close the figure."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / f"{stem}.pdf"
    fig.savefig(path, facecolor="white",
                metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    return path
