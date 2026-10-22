#!/usr/bin/env python3
"""Regenerate every PDF panel used by main.tex from the CSV tables in data/.

Usage:  python scripts/make_all_figures.py
Writes figures/*.pdf and prints ALL_FIGURES_REGENERATED_OK on success.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
SCRIPTS = (
    "fig01_pg_symbol_filters.py",
    "fig02_jacobian_h3_scaling.py",
    "fig03_controlled_damping_tradeoff.py",
    "fig04_highthroughput_summary.py",
    "fig05_absolute_spatial_convergence.py",
    "fig07_accuracy_cost_absolute_reference.py",
)


def main() -> None:
    env = dict(os.environ, MPLBACKEND="Agg")
    for script in SCRIPTS:
        print(f"[figure] {script}", flush=True)
        subprocess.run([sys.executable, str(HERE / script)], cwd=HERE, env=env, check=True)
    print("ALL_FIGURES_REGENERATED_OK")


if __name__ == "__main__":
    main()
