#!/usr/bin/env python3
"""Recompute the fixed-state Jacobian h-scaling experiment and fitted laws."""
from __future__ import annotations

import argparse
import csv
import json
import platform
from pathlib import Path

import numpy as np
import scipy

from jacobian_core import (
    AbsKPetrovGalerkin,
    amplification_row,
    compacton,
    ga_rho,
    largest_magnitude_eigenvalues,
    s1,
    s2,
    s3,
)


def write_csv(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def run(root: Path) -> None:
    generated = root / "generated"
    generated.mkdir(parents=True, exist_ok=True)

    length = 100.0
    p = 3.0
    c0 = 1.0
    c1, c2 = 2.0, -1.0
    mu = 1.0e-3
    dt = 0.025
    meshes = (250, 500, 1000, 2000, 4000)
    states = {
        "separated": (40.0, 60.0),
        "overlap": (48.0, 52.0),
    }

    records: list[dict] = []
    for state, (x1, x2) in states.items():
        for n in meshes:
            model = AbsKPetrovGalerkin(n, length, p, c0, mu)
            x = np.arange(n, dtype=float) * model.h
            u = compacton(x, c1, x1, p, length)
            u += compacton(x, c2, x2, p, length)
            values = largest_magnitude_eigenvalues(model, u, k=20)
            value = values[int(np.argmax(np.abs(values)))]
            z = dt * value
            records.append(
                {
                    "state": state,
                    "N": n,
                    "h": model.h,
                    "max_abs_u": float(np.max(np.abs(u))),
                    "lambda_real": float(np.real(value)),
                    "lambda_imag": float(np.imag(value)),
                    "abs_lambda_max": float(abs(value)),
                    "abs_z_max": float(abs(z)),
                    "M1_amp": float(abs(s1(z))),
                    "M2_amp": float(abs(s2(z))),
                    "M3_amp": float(abs(s3(z))),
                    "GA_rho0_amp": ga_rho(z, 0.0),
                    "GA_rho05_amp": ga_rho(z, 0.5),
                }
            )

    write_csv(generated / "jacobian_mesh_scaling.csv", records)

    fits: list[dict] = []
    quantities = ("abs_lambda_max", "M2_amp", "GA_rho0_amp")
    for state in states:
        state_rows = [row for row in records if row["state"] == state]
        for fit_name, selected in (
            ("all", state_rows),
            ("fine_h<=0.1", [row for row in state_rows if row["h"] <= 0.1]),
        ):
            h = np.array([row["h"] for row in selected], dtype=float)
            for quantity in quantities:
                values = np.array([row[quantity] for row in selected], dtype=float)
                slope, intercept = np.polyfit(np.log(h), np.log(values), 1)
                fits.append(
                    {
                        "state": state,
                        "fit_range": fit_name,
                        "quantity": quantity,
                        "slope_vs_h": float(slope),
                        "intercept": float(intercept),
                    }
                )
    write_csv(generated / "jacobian_mesh_scaling_fits.csv", fits)

    metadata = {
        "parameters": {
            "p": p,
            "c1": c1,
            "c2": c2,
            "c0": c0,
            "length": length,
            "mu": mu,
            "dt": dt,
            "meshes": list(meshes),
            "states": {key: list(value) for key, value in states.items()},
            "eigen_count": 20,
            "fine_fit_meshes": [1000, 2000, 4000],
        },
        "generalized_alpha_source": "ga_matrix/ga_rho in jacobian_core.py",
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "command": "python recompute_mesh_scaling.py",
    }
    (generated / "jacobian_mesh_scaling_config.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print("JACOBIAN_MESH_RECOMPUTE_OK")
    for state in states:
        fine = [
            row
            for row in fits
            if row["state"] == state and row["fit_range"] == "fine_h<=0.1"
        ]
        vals = {row["quantity"]: row["slope_vs_h"] for row in fine}
        print(
            f"{state}: lambda={vals['abs_lambda_max']:.12f} "
            f"M2={vals['M2_amp']:.12f} GA0={vals['GA_rho0_amp']:.12f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Jacobian reproducibility directory",
    )
    args = parser.parse_args()
    run(args.root.resolve())


if __name__ == "__main__":
    main()
