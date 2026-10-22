#!/usr/bin/env python3
"""Audit the full configuration -> states -> spectra -> filters -> Figure 2 chain."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from jacobian_core import ga_matrix


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_full(root: Path) -> None:
    env = os.environ.copy()
    for key in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        env[key] = "1"
    for script in ("recompute_collision_jacobian.py", "recompute_mesh_scaling.py"):
        subprocess.run(
            [sys.executable, str(root / script), "--root", str(root)],
            check=True,
            cwd=root,
            env=env,
        )


def compare_collision(root: Path) -> dict:
    frozen = read_csv(root / "frozen" / "jacobian_highfreq_modes_and_filters.csv")
    generated = read_csv(root / "generated" / "jacobian_highfreq_modes_and_filters.csv")
    if len(frozen) != 180 or len(generated) != 180:
        raise AssertionError(f"collision row count mismatch: {len(frozen)}, {len(generated)}")

    report: dict[str, dict] = {}
    for t in (5.0, 20.0, 30.0):
        frows = [row for row in frozen if float(row["t"]) == t]
        grows = [row for row in generated if float(row["t"]) == t]
        fvals = np.array(
            [complex(float(row["lambda_real"]), float(row["lambda_imag"])) for row in frows]
        )
        gvals = np.array(
            [complex(float(row["lambda_real"]), float(row["lambda_imag"])) for row in grows]
        )
        cost = np.abs(gvals[:, None] - fvals[None, :])
        gi, fi = linear_sum_assignment(cost)
        matched = cost[gi, fi]
        rel = matched / np.maximum(1.0, np.abs(fvals[fi]))
        amp_errors = {}
        for field in ("M1_amp", "M2_amp", "M3_amp", "GA_rho_0p5", "GA_rho_0"):
            gv = np.array([float(grows[i][field]) for i in gi])
            fv = np.array([float(frows[i][field]) for i in fi])
            amp_errors[field] = {
                "max_abs": float(np.max(np.abs(gv - fv))),
                "max_rel": float(np.max(np.abs(gv - fv) / np.maximum(1.0e-300, np.abs(fv)))),
            }
        max_abs = float(np.max(matched))
        max_rel = float(np.max(rel))
        if max_rel > 5.0e-11:
            raise AssertionError(f"collision spectrum mismatch at t={t}: {max_rel:.3e}")
        if max(err["max_rel"] for err in amp_errors.values()) > 5.0e-11:
            raise AssertionError(f"collision filter mismatch at t={t}")
        report[f"t={t:g}"] = {
            "max_complex_abs_error": max_abs,
            "max_complex_relative_error": max_rel,
            "amplification_errors": amp_errors,
        }
    return report


def keyed(rows: list[dict], fields: tuple[str, ...]) -> dict[tuple[str, ...], dict]:
    return {tuple(row[field] for field in fields): row for row in rows}


def compare_mesh(root: Path) -> dict:
    frozen = keyed(
        read_csv(root / "frozen" / "jacobian_mesh_scaling.csv"),
        ("state", "N"),
    )
    generated = keyed(
        read_csv(root / "generated" / "jacobian_mesh_scaling.csv"),
        ("state", "N"),
    )
    if frozen.keys() != generated.keys():
        raise AssertionError("mesh keys differ")

    fields = (
        "max_abs_u",
        "abs_lambda_max",
        "abs_z_max",
        "M1_amp",
        "M2_amp",
        "M3_amp",
        "GA_rho0_amp",
        "GA_rho05_amp",
    )
    field_errors = {}
    for field in fields:
        abs_errors, rel_errors = [], []
        for key in frozen:
            a = float(frozen[key][field])
            b = float(generated[key][field])
            abs_errors.append(abs(a - b))
            rel_errors.append(abs(a - b) / max(1.0e-300, abs(a)))
        field_errors[field] = {
            "max_abs": max(abs_errors),
            "max_rel": max(rel_errors),
        }
    # ARPACK may choose either member of a conjugate pair; compare magnitudes.
    for field in ("lambda_real", "lambda_imag"):
        errors = [
            abs(abs(float(frozen[key][field])) - abs(float(generated[key][field])))
            for key in frozen
        ]
        field_errors[f"abs({field})"] = {"max_abs": max(errors)}

    if max(field_errors[field]["max_rel"] for field in fields) > 5.0e-11:
        raise AssertionError("mesh data exceed relative tolerance")

    ffits = keyed(
        read_csv(root / "frozen" / "jacobian_mesh_scaling_fits.csv"),
        ("state", "fit_range", "quantity"),
    )
    gfits = keyed(
        read_csv(root / "generated" / "jacobian_mesh_scaling_fits.csv"),
        ("state", "fit_range", "quantity"),
    )
    fit_errors = {}
    for field in ("slope_vs_h", "intercept"):
        error = max(abs(float(ffits[key][field]) - float(gfits[key][field])) for key in ffits)
        fit_errors[field] = error
        if error > 1.0e-10:
            raise AssertionError(f"fit mismatch in {field}: {error:.3e}")

    return {"data_errors": field_errors, "fit_errors": fit_errors}


def compare_ga_matrix(root: Path) -> dict:
    ga_path = root.parent / "generalized_alpha" / "ga_core.py"
    spec = importlib.util.spec_from_file_location("corrected_ga_core", ga_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {ga_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    maximum = 0.0
    for rho in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0):
        for z in (-0.2, 0.1 + 3.0j, -4.0 + 800.0j, -0.3 - 1.0e5j):
            maximum = max(
                maximum,
                float(np.max(np.abs(ga_matrix(z, rho) - module.ga_matrix(z, rho)))),
            )
    if maximum > 2.0e-15:
        raise AssertionError(f"ga_matrix copies differ: {maximum:.3e}")
    return {"max_matrix_difference": maximum, "corrected_source": str(ga_path)}


def verify_figure_sources(root: Path) -> dict:
    figure_root = root.parents[1] / "figure_sources"
    fig1 = (figure_root / "fig01_pg_symbol_filters.py").read_text(encoding="utf-8")
    fig2 = (figure_root / "fig02_jacobian_h3_scaling.py").read_text(encoding="utf-8")
    combined = fig1 + "\n" + fig2
    if "step_GA" in combined or "step_GA_fast" in combined:
        raise AssertionError("Figure 1/2 source calls a nonlinear GA step")
    if "from jacobian_core import ga_rho" not in fig1:
        raise AssertionError("Figure 1 does not import the canonical ga_rho")
    if "midpoint_amp" in fig2 or "'M1_amp'" not in fig2:
        raise AssertionError("Figure 2 does not use the regenerated M1_amp column")

    figure_data = read_csv(figure_root / "data" / "fig2_jacobian_h3_scaling.csv")
    generated = read_csv(root / "generated" / "jacobian_mesh_scaling.csv")
    fkey = keyed(figure_data, ("state", "N"))
    gkey = keyed(generated, ("state", "N"))
    max_rel = 0.0
    for key in fkey:
        for field in (
            "abs_lambda_max",
            "M1_amp",
            "M2_amp",
            "M3_amp",
            "GA_rho0_amp",
            "GA_rho05_amp",
        ):
            a, b = float(fkey[key][field]), float(gkey[key][field])
            max_rel = max(max_rel, abs(a - b) / max(1.0e-300, abs(a)))
    if max_rel > 5.0e-11:
        raise AssertionError(f"Figure 2 data not traced to generated mesh data: {max_rel:.3e}")
    return {
        "figure_1_generalized_alpha": "canonical ga_matrix/ga_rho only",
        "figure_2_max_relative_data_error": max_rel,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parent
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="rerun the M2 trajectory, both ARPACK experiments, and then audit",
    )
    parser.add_argument(
        "--_in-place-full", action="store_true", help=argparse.SUPPRESS
    )
    args = parser.parse_args()
    root = args.root.resolve()

    # ARPACK may return machine-dependent last-bit variations.  A full audit
    # therefore runs in a disposable package copy and never rewrites the frozen
    # generated tables in the user's working tree.
    if args.full and not args._in_place_full:
        package_root = root.parents[1]
        with tempfile.TemporaryDirectory(prefix="compacton-jacobian-full-") as tmp:
            sandbox_package = Path(tmp) / package_root.name
            shutil.copytree(package_root, sandbox_package)
            sandbox_root = sandbox_package / "reproducibility" / "jacobian"
            env = dict(os.environ)
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            for key in (
                "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
            ):
                env[key] = "1"
            subprocess.check_call(
                [
                    sys.executable,
                    str(sandbox_root / "verify_jacobian_traceability.py"),
                    "--root",
                    str(sandbox_root),
                    "--full",
                    "--_in-place-full",
                ],
                cwd=sandbox_root,
                env=env,
            )
        print("JACOBIAN_FULL_AUDIT_DISPOSABLE_COPY_OK")
        return

    if args.full:
        run_full(root)

    report = {
        "status": "JACOBIAN_TRACEABILITY_OK",
        "provenance": {
            "driver_status": "reconstructed authoritative drivers",
            "original_historical_driver_recovered_byte_for_byte": False,
            "basis": "archived equations, parameters, reports, frozen CSV outputs, and independent full recomputation",
        },
        "collision": compare_collision(root),
        "mesh_scaling": compare_mesh(root),
        "generalized_alpha_matrix": compare_ga_matrix(root),
        "figure_sources": verify_figure_sources(root),
        "artifacts": {},
    }
    for relative in (
        "raw/jacobian_collision_snapshots.npz",
        "generated/jacobian_highfreq_modes_and_filters.csv",
        "generated/jacobian_mesh_scaling.csv",
        "generated/jacobian_mesh_scaling_fits.csv",
        "frozen/jacobian_highfreq_modes_and_filters.csv",
        "frozen/jacobian_mesh_scaling.csv",
        "frozen/jacobian_mesh_scaling_fits.csv",
    ):
        path = root / relative
        report["artifacts"][relative] = {"sha256": sha256(path), "bytes": path.stat().st_size}

    json_path = root / "generated" / "JACOBIAN_VERIFICATION.json"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    text = [
        "JACOBIAN_TRACEABILITY_OK",
        "",
        "Provenance:",
        "  authoritative drivers were reconstructed from the documented equations and parameters",
        "  the original historical eigenvalue drivers were not recovered byte for byte",
        "  full recomputation matches every frozen numerical output within the stated tolerances",
        "",
        "Collision trajectory and spectra:",
    ]
    for time_key, values in report["collision"].items():
        text.append(
            f"  {time_key}: max relative eigenvalue error = "
            f"{values['max_complex_relative_error']:.3e}"
        )
    text.extend(
        [
            "",
            "Fixed-state mesh scaling:",
            f"  max slope error = {report['mesh_scaling']['fit_errors']['slope_vs_h']:.3e}",
            f"  max intercept error = {report['mesh_scaling']['fit_errors']['intercept']:.3e}",
            "",
            "Generalized-alpha matrix:",
            f"  max difference between canonical and corrected campaign copies = "
            f"{report['generalized_alpha_matrix']['max_matrix_difference']:.3e}",
            "",
            "Figure sources:",
            "  Figure 1 computes generalized-alpha amplification only through ga_matrix/ga_rho.",
            "  Figure 2 reads data regenerated by recompute_mesh_scaling.py.",
            "  No Figure 1/2 source calls step_GA or step_GA_fast.",
            "",
        ]
    )
    txt_path = root / "generated" / "JACOBIAN_VERIFICATION.txt"
    txt_path.write_text("\n".join(text), encoding="utf-8")
    print(txt_path.read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
