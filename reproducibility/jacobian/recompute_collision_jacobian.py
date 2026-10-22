#!/usr/bin/env python3
"""Recompute the nonlinear-collision Jacobian experiment from initial data.

Outputs
-------
raw/jacobian_collision_snapshots.npz
    x grid and U(t) at t=5,20,30.
generated/jacobian_highfreq_modes_and_filters.csv
    60 largest-modulus eigenvalues at every snapshot and their time-filter
    amplification factors.
generated/jacobian_collision_summary.csv
    One-row-per-snapshot summary used by the manuscript audit.
generated/jacobian_collision_config.json
    Complete parameter and software metadata.
"""
from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import scipy

from jacobian_core import (
    AbsKPetrovGalerkin,
    CollisionConfig,
    amplification_row,
    compacton,
    largest_magnitude_eigenvalues,
    step_m2,
)


def write_csv(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def run(root: Path) -> None:
    config = CollisionConfig()
    generated = root / "generated"
    raw = root / "raw"
    generated.mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)

    model = AbsKPetrovGalerkin(
        config.n, config.length, config.p, config.c0, config.mu
    )
    x = np.arange(config.n, dtype=float) * model.h
    u = compacton(x, config.c1, config.x1, config.p, config.length)
    u += compacton(x, config.c2, config.x2, config.p, config.length)
    u_previous = None

    snapshot_steps = {
        int(round(t / config.dt)): float(t) for t in config.snapshot_times
    }
    if any(abs(step * config.dt - t) > 1.0e-13 for step, t in snapshot_steps.items()):
        raise RuntimeError("snapshot times must be exact time-step multiples")

    snapshots: dict[float, np.ndarray] = {}
    step_records: list[dict] = []
    start = time.perf_counter()
    final_step = max(snapshot_steps)
    for step in range(1, final_step + 1):
        x_new, iterations, ok, residual = step_m2(
            model,
            u,
            config.dt,
            u_previous,
            tol=config.newton_tol,
            maxit=config.newton_maxit,
        )
        if not ok:
            raise RuntimeError(
                f"M2 Newton failure at step={step}, t={step*config.dt:.12g}, "
                f"iterations={iterations}, residual={residual:.6e}"
            )
        u_previous, u = u, x_new
        if step in snapshot_steps:
            t = snapshot_steps[step]
            snapshots[t] = u.copy()
            step_records.append(
                {
                    "t": t,
                    "step": step,
                    "max_abs_u": float(np.max(np.abs(u))),
                    "newton_iterations": int(iterations),
                    "newton_residual_inf": float(residual),
                }
            )

    snapshot_payload = {"x": x}
    for t in config.snapshot_times:
        snapshot_payload[f"U_t{t:g}"] = snapshots[t]
    np.savez(raw / "jacobian_collision_snapshots.npz", **snapshot_payload)

    rows: list[dict] = []
    summary: list[dict] = []
    stage_name = {5.0: "pre-collision", 20.0: "collision", 30.0: "post-collision"}
    for t in config.snapshot_times:
        values = largest_magnitude_eigenvalues(
            model, snapshots[t], config.eigen_count
        )
        amp_rows = []
        for mode_index, value in enumerate(values):
            amps = amplification_row(complex(value), config.dt)
            row = {
                "t": t,
                "mode_index": mode_index,
                "lambda_real": float(np.real(value)),
                "lambda_imag": float(np.imag(value)),
                **amps,
            }
            rows.append(row)
            amp_rows.append(row)

        summary.append(
            {
                "t": t,
                "stage": stage_name[t],
                "max_abs_u": float(np.max(np.abs(snapshots[t]))),
                "max_abs_Im_dt_lambda": max(abs(r["z_imag"]) for r in amp_rows),
                "median_M1_amp": float(np.median([r["M1_amp"] for r in amp_rows])),
                "median_M2_amp": float(np.median([r["M2_amp"] for r in amp_rows])),
                "median_M3_amp": float(np.median([r["M3_amp"] for r in amp_rows])),
                "median_GA_rho_0p5": float(
                    np.median([r["GA_rho_0p5"] for r in amp_rows])
                ),
                "median_GA_rho_0": float(
                    np.median([r["GA_rho_0"] for r in amp_rows])
                ),
            }
        )

    write_csv(generated / "jacobian_highfreq_modes_and_filters.csv", rows)
    write_csv(generated / "jacobian_collision_summary.csv", summary)
    write_csv(generated / "jacobian_collision_step_diagnostics.csv", step_records)

    metadata = {
        "configuration": asdict(config),
        "spatial_formula": (
            "J(U)=-A^{-1}[mu*h^2*D-c0*B+(B+C)diag(p*abs(U)^(p-1))]"
        ),
        "trajectory_method": "Padé R_{1,2} fractional step (M2)",
        "eigensolver": {
            "library": "scipy.sparse.linalg.eigs (ARPACK)",
            "which": "LM",
            "tolerance": 1.0e-10,
            "maxiter": 50000,
            "starting_vector": "deterministic_v0 in jacobian_core.py",
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "elapsed_seconds": time.perf_counter() - start,
        "command": "python recompute_collision_jacobian.py",
    }
    (generated / "jacobian_collision_config.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print("JACOBIAN_COLLISION_RECOMPUTE_OK")
    for row in summary:
        print(
            f"t={row['t']:g} max_abs_u={row['max_abs_u']:.12g} "
            f"max_abs_Im_dt_lambda={row['max_abs_Im_dt_lambda']:.12g} "
            f"median_M2={row['median_M2_amp']:.12g}"
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
