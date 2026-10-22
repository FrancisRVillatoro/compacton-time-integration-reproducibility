"""Common CSV reader and consistent single-panel geometry (72 x 54 mm axes)."""
from __future__ import annotations
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from figure_style_final import MM, AX_W, AX_H, H, BOTTOM, style, save

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
OUT = ROOT / 'figures'


def rows(name):
    with (DATA / name).open('r', newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def column(records, field):
    return np.array([float(r[field]) for r in records], dtype=float)


def single_panel():
    # The current manuscript embeds these panels separately at 49% line width.
    # Use the same 72x54 mm physical plotting box in every panel.
    width_mm = 96.0
    fig = plt.figure(figsize=(width_mm * MM, H))
    ax = fig.add_axes([18.0 / width_mm, BOTTOM / H,
                       AX_W / (width_mm * MM), AX_H / H])
    return fig, style(ax)


def output(fig, name):
    save(fig, OUT / name)
