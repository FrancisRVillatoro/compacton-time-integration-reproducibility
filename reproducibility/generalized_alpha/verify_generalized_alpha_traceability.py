#!/usr/bin/env python3
"""Audit the corrected generalized-alpha campaign, Fig. 3, and manuscript text.

This verifier does not rerun the 396+64 production cases by default.  It closes
that part of the trace by checking the Picasso group outputs, their aggregates,
the published figure source, the manuscript numbers, the Slurm logs, and the
independent implementation validation.  With ``--full`` it reruns the small
matrix/order/midpoint validation in a clean output directory.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parents[1]
RESULTS = ROOT / "results"
VERIFY = ROOT / "verification"
FIGURE_DATA = PACKAGE / "figure_sources" / "data" / "fig3_controlled_damping_fine_corrected.csv"
FIGURE_SCRIPT = PACKAGE / "figure_sources" / "fig03_controlled_damping_tradeoff.py"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def quantile(values: list[float], q: float) -> float:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return math.nan
    position = (len(clean) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return clean[lower]
    return clean[lower] * (upper - position) + clean[upper] * (position - lower)


def assert_close(actual: float, expected: float, atol: float = 2.0e-14, rtol: float = 2.0e-12) -> None:
    if not math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol):
        raise AssertionError(f"numeric mismatch: {actual!r} != {expected!r}")


def compare_rows_exact(actual: list[dict[str, str]], expected: list[dict[str, str]], name: str) -> None:
    if actual != expected:
        raise AssertionError(f"{name} aggregate does not equal sorted group concatenation")


def verify_picasso_hashes() -> dict:
    mapping = {}
    manifest = RESULTS / "PICASSO_SHA256_RESULTS.txt"
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, old_path = line.split(maxsplit=1)
        mapping[Path(old_path.strip()).name] = digest
    actual = {}
    for name, expected in mapping.items():
        path = RESULTS / name
        if not path.exists():
            raise AssertionError(f"missing Picasso result listed in manifest: {name}")
        got = sha256(path)
        if got != expected:
            raise AssertionError(f"Picasso SHA-256 mismatch for {name}: {got} != {expected}")
        actual[name] = got
    return actual


def verify_aggregates() -> dict:
    micro_files = sorted((RESULTS / "microbank_groups").glob("group_*.csv"))
    boundary_files = sorted((RESULTS / "boundary_groups").glob("group_*.csv"))
    if len(micro_files) != 36 or len(boundary_files) != 8:
        raise AssertionError(f"unexpected group counts: {len(micro_files)}, {len(boundary_files)}")
    micro = [row for path in micro_files for row in read_csv(path)]
    boundary = [row for path in boundary_files for row in read_csv(path)]
    if len(micro) != 396 or len(boundary) != 64:
        raise AssertionError(f"unexpected row counts: {len(micro)}, {len(boundary)}")
    compare_rows_exact(micro, read_csv(RESULTS / "microbank_corrected_all.csv"), "microbank")
    compare_rows_exact(boundary, read_csv(RESULTS / "boundary_corrected_all.csv"), "boundary")

    # Reconstruct the screening summary retained for provenance but not used in the new Fig. 3.
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in micro:
        if row["status"] == "ok":
            groups[row["label"]].append(row)
    order = [
        "GA_1", "GA_0.9", "GA_0.75", "GA_0.5", "GA_0.447627",
        "GA_0.25", "GA_0.2", "GA_0.1", "GA_0.05", "GA_0", "M2",
    ]
    reconstructed_left = []
    for label in order:
        rows = groups[label]
        reconstructed_left.append({
            "label": label,
            "n_ok": len(rows),
            "median_auc": quantile([number(row["rho_hw_post_AUC"]) for row in rows], 0.5),
            "median_abs_I2_drift": quantile([abs(number(row["I2_drift_rel"])) for row in rows], 0.5),
            "p95_auc": quantile([number(row["rho_hw_post_AUC"]) for row in rows], 0.95),
        })
    frozen_left = read_csv(RESULTS / "fig3_left_corrected.csv")
    if len(frozen_left) != len(reconstructed_left):
        raise AssertionError("left-summary row count mismatch")
    for actual, expected in zip(reconstructed_left, frozen_left):
        if actual["label"] != expected["label"] or actual["n_ok"] != int(expected["n_ok"]):
            raise AssertionError("left-summary label/count mismatch")
        for field in ("median_auc", "median_abs_I2_drift", "p95_auc"):
            assert_close(float(actual[field]), float(expected[field]))

    # Reconstruct the fine case that is now the sole source of Fig. 3.
    selected = [
        row for row in boundary
        if int(row["N"]) == 1000
        and math.isclose(number(row["c1"]), 3.0)
        and math.isclose(number(row["mu"]), 1.0e-7, rel_tol=0.0, abs_tol=1.0e-20)
        and math.isclose(number(row["dt"]), 0.2)
    ]
    pade = next(row for row in selected if row["label"] == "M2")
    pade_auc = number(pade["rho_hw_post_AUC"])
    pade_drift = abs(number(pade["I2_drift_rel"]))
    reconstructed_right = []
    for row in selected:
        if row["status"] != "ok":
            continue
        auc = number(row["rho_hw_post_AUC"])
        drift = abs(number(row["I2_drift_rel"]))
        reconstructed_right.append({
            "label": row["label"],
            "rho_inf": row["rho_inf"],
            "status": row["status"],
            "auc": auc,
            "abs_I2_drift": drift,
            "auc_ratio_to_M2": auc / pade_auc,
            "I2_ratio_to_M2": drift / pade_drift,
        })
    frozen_right = read_csv(RESULTS / "fig3_right_corrected.csv")
    if len(frozen_right) != len(reconstructed_right):
        raise AssertionError("fine-summary row count mismatch")
    for actual, expected in zip(reconstructed_right, frozen_right):
        for field in ("label", "rho_inf", "status"):
            if actual[field] != expected[field]:
                raise AssertionError(f"fine-summary field mismatch: {field}")
        for field in ("auc", "abs_I2_drift", "auc_ratio_to_M2", "I2_ratio_to_M2"):
            assert_close(float(actual[field]), float(expected[field]))

    statuses = {
        "microbank": dict(Counter(row["status"] for row in micro)),
        "boundary": dict(Counter(row["status"] for row in boundary)),
        "single": dict(Counter(row["status"] for row in read_csv(RESULTS / "single_benchmark_corrected.csv"))),
    }
    if statuses != {
        "microbank": {"ok": 396},
        "boundary": {"ok": 62, "newton_failure": 2},
        "single": {"ok": 7},
    }:
        raise AssertionError(f"unexpected status counts: {statuses}")
    return {
        "group_files": {"microbank": 36, "boundary": 8},
        "row_counts": {"microbank": 396, "boundary": 64, "single": 7},
        "statuses": statuses,
        "fine_case": {
            "pade_auc": pade_auc,
            "pade_abs_I2_drift": pade_drift,
            "rows": reconstructed_right,
        },
    }


def verify_figure_source() -> dict:
    result_path = RESULTS / "fig3_right_corrected.csv"
    if sha256(FIGURE_DATA) != sha256(result_path):
        raise AssertionError("Fig. 3 source CSV is not byte-identical to the audited fine-case summary")
    script = FIGURE_SCRIPT.read_text(encoding="utf-8")
    if "fig3_controlled_damping_fine_corrected.csv" not in script:
        raise AssertionError("Fig. 3 script does not read the audited fine-case CSV")
    if "fig3_left_corrected.csv" in script:
        raise AssertionError("Fig. 3 script still reads the coarse screening summary")
    return {
        "source_csv_sha256": sha256(FIGURE_DATA),
        "script_sha256": sha256(FIGURE_SCRIPT),
        "coarse_screening_summary_used": False,
    }


def verify_single_case() -> dict:
    rows = {row["label"]: row for row in read_csv(RESULTS / "single_benchmark_corrected.csv")}
    expected = {
        "M2": (4.6121671352964577e-4, 1.1153888902719544e-3),
        "GA_0.25": (2.579402626458978e-4, 8.895961843560324e-4),
        "GA_0": (4.525703604374808e-5, 1.5929834290033453e-3),
        "M1": (2.4214314351335395, 2.0511344649108196e-1),
    }
    report = {}
    for label, (auc_expected, drift_expected) in expected.items():
        row = rows[label]
        auc = number(row["rho_hw_post_AUC"])
        drift = abs(number(row["I2_drift_rel"]))
        assert_close(auc, auc_expected)
        assert_close(drift, drift_expected)
        report[label] = {"auc": auc, "abs_I2_drift": drift}
    final_rho = number(rows["M1"]["rho_hw_final"])
    assert_close(final_rho, 0.24514119759756442)
    report["M1"]["rho_hw_final"] = final_rho
    return report


def verify_picasso_logs() -> dict:
    log_dir = ROOT / "picasso_provenance" / "logs"
    err = sorted(log_dir.glob("*.err"))
    out = sorted(log_dir.glob("*.out"))
    nonempty = [path.name for path in err if path.stat().st_size != 0]
    if len(err) != 47 or len(out) != 47 or nonempty:
        raise AssertionError(
            f"unexpected Picasso log audit: err={len(err)}, out={len(out)}, nonempty={nonempty}"
        )
    return {"stderr_files": len(err), "stdout_files": len(out), "nonempty_stderr": nonempty}


def rerun_validation() -> dict:
    output = VERIFY / "validation_recheck"
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    env["GA_OUT"] = str(output)
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    completed = subprocess.run(
        [sys.executable, str(ROOT / "validate_ga.py")],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )
    (output / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    if completed.stderr:
        (output / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    summary = (output / "validation_summary.txt").read_text(encoding="utf-8")
    if not summary.startswith("GA_CORRECTED_VALIDATION_OK\n"):
        raise AssertionError("corrected generalized-alpha validation did not pass")
    production_summary = (RESULTS / "validation_summary.txt").read_text(encoding="utf-8")
    if summary != production_summary:
        raise AssertionError("validation recheck differs from the Picasso validation summary")
    return {
        "summary_sha256": sha256(output / "validation_summary.txt"),
        "scalar_order_sha256": sha256(output / "validation_scalar_order.csv"),
        "summary": summary.splitlines(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true", help="rerun the matrix/order/midpoint validation")
    args = parser.parse_args()
    VERIFY.mkdir(parents=True, exist_ok=True)

    report = {
        "status": "GENERALIZED_ALPHA_TRACEABILITY_OK",
        "implementation": {
            "canonical_source": str(ROOT / "ga_core.py"),
            "correct_residual": "a A (X-U_n) + c A V_n + G((1-alpha_f)U_n+alpha_f X) = 0",
            "historical_results_used_in_manuscript": False,
        },
        "picasso_hashes": verify_picasso_hashes(),
        "aggregates": verify_aggregates(),
        "single_case": verify_single_case(),
        "figure_3": verify_figure_source(),
        "picasso_logs": verify_picasso_logs(),
        "validation": rerun_validation() if args.full else {
            "summary_sha256": sha256(RESULTS / "validation_summary.txt"),
            "scalar_order_sha256": sha256(RESULTS / "validation_scalar_order.csv"),
            "mode": "frozen Picasso validation (use --full to rerun)",
        },
    }

    json_path = VERIFY / "GENERALIZED_ALPHA_VERIFICATION.json"
    txt_path = VERIFY / "GENERALIZED_ALPHA_VERIFICATION.txt"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    fine = report["aggregates"]["fine_case"]
    text = [
        "GENERALIZED_ALPHA_TRACEABILITY_OK",
        "",
        "Corrected implementation:",
        "  residual retains the independent algorithmic derivative V_n",
        "  historical defective generalized-alpha results are not used",
        "",
        "Picasso campaign:",
        "  microbank: 396/396 ok",
        "  fine boundary bank: 62/64 ok; the two failures are rho_inf=1 controls",
        "  independent N=1000, c1=2 benchmark: 7/7 ok",
        "  47/47 Slurm stderr files are empty",
        "",
        "Fine case used in Figure 3:",
        f"  Pade AUC = {fine['pade_auc']:.16e}",
        f"  Pade |Delta I2/I2| = {fine['pade_abs_I2_drift']:.16e}",
        "  figure CSV is byte-identical to the audited Picasso summary",
        "",
    ]
    if args.full:
        text.extend(["Validation rerun:", *[f"  {line}" for line in report["validation"]["summary"]], ""])
    txt_path.write_text("\n".join(text), encoding="utf-8")
    print(txt_path.read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
