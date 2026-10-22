#!/usr/bin/env python3
"""Audit the recovered Banks 01--04 configuration and production chain.

Run from the package root:

    python scripts/audit_robustness_banks.py

The script reads the complete recovered case and summary CSVs, checks the
MATLAB/SLURM provenance, and rewrites the two compact audit tables in data/.
It intentionally uses only the Python standard library.
"""

from __future__ import annotations

import csv
import math
from collections import Counter
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
PROV = ROOT / "provenance" / "robustness_banks"
DATA = ROOT / "data"

BANKS = [
    {
        "name": "bank01_mu_dt",
        "label": "Primary",
        "cases": "bank01_mu_dt_cases.csv",
        "summary": "bank01_mu_dt_summary.csv",
        "expected": 1440,
        "post_rule": "60 (fixed observation-window start)",
    },
    {
        "name": "bank02_dx_refinement",
        "label": "Spatial refinement",
        "cases": "bank02_dx_refinement_cases.csv",
        "summary": "bank02_dx_refinement_summary.csv",
        "expected": 1296,
        "post_rule": "60 (fixed observation-window start)",
    },
    {
        "name": "bank03_collision_severity",
        "label": "Collision severity",
        "cases": "bank03_collision_severity_cases.csv",
        "summary": "bank03_collision_severity_summary.csv",
        "expected": 1728,
        "post_rule": "min(0.8*Tmax, max(10, (x02-x01)/abs(c1-c2) + 10))",
    },
    {
        "name": "bank04_long_time",
        "label": "Long time",
        "cases": "bank04_long_time_cases.csv",
        "summary": "bank04_long_time_summary.csv",
        "expected": 576,
        "post_rule": "60 (fixed observation-window start)",
    },
]

SETTINGS_FIELDS = [
    "bank",
    "label",
    "n_cases",
    "method",
    "scenario",
    "p",
    "dx",
    "dt",
    "mu",
    "Lmax",
    "Tmax",
    "c0",
    "c1",
    "c2",
    "x01",
    "x02",
    "post_t0",
    "save_raw",
    "save_snapshots",
    "nsave_target",
    "newton_tol",
    "newton_max_iter",
    "blowup_factor",
    "post_t0_rule",
    "newton_stopping_test",
    "newton_cap_semantics",
    "amplitude_stop_condition",
]

CAP_FIELDS = [
    "bank",
    "method",
    "assigned",
    "status_ok",
    "status_ok_with_any_cap_hit",
    "max_total_cap_hits_in_one_status_ok_run",
    "max_stored_newton_updates_among_status_ok",
    "note",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def numeric_key(value: str) -> tuple[int, float | str]:
    try:
        return (0, float(value))
    except ValueError:
        return (1, value)


def compact_unique(rows: Iterable[dict[str, str]], field: str) -> str:
    values = sorted({row[field] for row in rows}, key=numeric_key)
    return ";".join(values)


def as_float(value: str) -> float:
    return float(value.strip())


def assert_all_numeric(rows: list[dict[str, str]], field: str, expected: float) -> None:
    bad = [row.get("run_id", "?") for row in rows if not math.isclose(as_float(row[field]), expected, rel_tol=0.0, abs_tol=1e-15)]
    if bad:
        raise AssertionError(f"{field} != {expected} in run ids {bad[:10]}")


def audit_sources() -> None:
    matlab = PROV / "matlab"
    slurm = PROV / "slurm"

    dispatcher = (matlab / "run_absK_case.m").read_text(encoding="utf-8")
    required_dispatch = [
        "readtable(cases_csv",
        "meta.post_t0 = row.post_t0;",
        "meta.newton_tol = row.newton_tol;",
        "meta.newton_max_iter = row.newton_max_iter;",
        "meta.blowup_factor = row.blowup_factor;",
        "deFrutosAbsDissipative_HPC",
        "Pade21enT_HPC",
        "kpqSimulacionMP4T_HPC",
    ]
    for marker in required_dispatch:
        if marker not in dispatcher:
            raise AssertionError(f"dispatcher marker missing: {marker}")

    drivers = [
        matlab / "deFrutosAbsDissipative_HPC.m",
        matlab / "Pade21enT_HPC.m",
        matlab / "kpqSimulacionMP4T_HPC.m",
    ]
    for path in drivers:
        text = path.read_text(encoding="utf-8")
        required = [
            "TOL = meta.newton_tol;",
            "newtIter < meta.newton_max_iter",
            "newtonStats.fail_count = newtonStats.fail_count + 1;",
            "max(abs(Un)) > meta.blowup_factor*meta.max_abs_u0",
        ]
        for marker in required:
            if marker not in text:
                raise AssertionError(f"{path.name}: marker missing: {marker}")
        if "norm(Uk-Ukp1) > TOL" not in text and "norm(UVk-UVkp1) > TOL" not in text:
            raise AssertionError(f"{path.name}: successive-iterate Newton test not found")

    production_launchers = [
        slurm / "run_bank01_mu_dt.slurm",
        slurm / "run_bank02_dx_refinement.slurm",
        slurm / "run_bank03_collision_severity.slurm",
        slurm / "run_bank04_long_time.slurm",
    ]
    for path in production_launchers:
        text = path.read_text(encoding="utf-8")
        for marker in ("module load matlab/R2024b", "run_absK_case", "matlab -singleCompThread"):
            if marker not in text:
                raise AssertionError(f"{path.name}: launch marker missing: {marker}")


def main() -> None:
    audit_sources()
    DATA.mkdir(parents=True, exist_ok=True)

    setting_rows: list[dict[str, str | int]] = []
    cap_rows: list[dict[str, str | int]] = []

    for bank in BANKS:
        cases_path = PROV / "cases" / bank["cases"]
        summary_path = PROV / "summaries" / bank["summary"]
        cases = read_csv(cases_path)
        summary = read_csv(summary_path)

        if len(cases) != bank["expected"] or len(summary) != bank["expected"]:
            raise AssertionError(
                f"{bank['name']}: rows cases/summary = {len(cases)}/{len(summary)}, expected {bank['expected']}"
            )

        case_ids = {row["run_id"] for row in cases}
        summary_ids = {row["run_id"] for row in summary}
        if case_ids != summary_ids:
            raise AssertionError(f"{bank['name']}: case/summary run IDs differ")

        assert_all_numeric(cases, "newton_tol", 1e-10)
        assert_all_numeric(cases, "newton_max_iter", 10.0)
        assert_all_numeric(cases, "blowup_factor", 50.0)

        if bank["name"] in {"bank01_mu_dt", "bank02_dx_refinement", "bank04_long_time"}:
            assert_all_numeric(cases, "post_t0", 60.0)
        else:
            for row in cases:
                tmax = as_float(row["Tmax"])
                x01 = as_float(row["x01"])
                x02 = as_float(row["x02"])
                c1 = as_float(row["c1"])
                c2 = as_float(row["c2"])
                expected_post = min(0.8 * tmax, max(10.0, (x02 - x01) / abs(c1 - c2) + 10.0))
                if not math.isclose(as_float(row["post_t0"]), expected_post, rel_tol=0.0, abs_tol=1e-12):
                    raise AssertionError(
                        f"{bank['name']} run {row['run_id']}: post_t0 {row['post_t0']} != {expected_post}"
                    )

        setting_rows.append(
            {
                "bank": bank["name"],
                "label": bank["label"],
                "n_cases": len(cases),
                "method": compact_unique(cases, "method"),
                "scenario": compact_unique(cases, "scenario"),
                "p": compact_unique(cases, "p"),
                "dx": compact_unique(cases, "dx"),
                "dt": compact_unique(cases, "dt"),
                "mu": compact_unique(cases, "mu"),
                "Lmax": compact_unique(cases, "Lmax"),
                "Tmax": compact_unique(cases, "Tmax"),
                "c0": compact_unique(cases, "c0"),
                "c1": compact_unique(cases, "c1"),
                "c2": compact_unique(cases, "c2"),
                "x01": compact_unique(cases, "x01"),
                "x02": compact_unique(cases, "x02"),
                "post_t0": compact_unique(cases, "post_t0"),
                "save_raw": compact_unique(cases, "save_raw"),
                "save_snapshots": compact_unique(cases, "save_snapshots"),
                "nsave_target": compact_unique(cases, "nsave_target"),
                "newton_tol": compact_unique(cases, "newton_tol"),
                "newton_max_iter": compact_unique(cases, "newton_max_iter"),
                "blowup_factor": compact_unique(cases, "blowup_factor"),
                "post_t0_rule": bank["post_rule"],
                "newton_stopping_test": "norm(iterate_new-iterate_old,2) <= newton_tol",
                "newton_cap_semantics": "newton_max_iter=10 is the loop cap; stored newton_max <= 9",
                "amplitude_stop_condition": "max(abs(U_n)) > blowup_factor * max(abs(U_0))",
            }
        )

        methods = sorted({row["method"] for row in cases})
        for method in methods:
            assigned = sum(row["method"] == method for row in cases)
            ok_rows = [row for row in summary if row["method"] == method and row["status"] == "ok"]
            with_cap = [row for row in ok_rows if int(float(row["newton_fail_count"])) > 0]
            cap_rows.append(
                {
                    "bank": bank["name"],
                    "method": method,
                    "assigned": assigned,
                    "status_ok": len(ok_rows),
                    "status_ok_with_any_cap_hit": len(with_cap),
                    "max_total_cap_hits_in_one_status_ok_run": max(
                        (int(float(row["newton_fail_count"])) for row in ok_rows), default=0
                    ),
                    "max_stored_newton_updates_among_status_ok": max(
                        (int(float(row["newton_max"])) for row in ok_rows), default=0
                    ),
                    "note": "cap hits are recorded separately and do not change status",
                }
            )

        status_counts = Counter(row["status"] for row in summary)
        print(f"{bank['name']}: cases={len(cases)}, summary={len(summary)}, status={dict(status_counts)}")

    settings_out = DATA / "robustness_bank_settings_recovered.csv"
    with settings_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SETTINGS_FIELDS)
        writer.writeheader()
        writer.writerows(setting_rows)

    cap_out = DATA / "robustness_bank_newton_cap_audit.csv"
    with cap_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CAP_FIELDS)
        writer.writeheader()
        writer.writerows(cap_rows)

    print(f"wrote {settings_out.relative_to(ROOT)}")
    print(f"wrote {cap_out.relative_to(ROOT)}")
    print("ROBUSTNESS_BANK_TRACEABILITY_OK")


if __name__ == "__main__":
    main()
