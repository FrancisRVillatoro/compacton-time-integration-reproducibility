#!/usr/bin/env python3
"""Cross-check v1.1.0 metadata and compact manuscript-facing tables."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    require((ROOT / "VERSION").read_text().strip() == "1.1.0", "VERSION is not 1.1.0")
    require("version: 1.1.0" in (ROOT / "CITATION.cff").read_text(), "CITATION.cff version mismatch")
    require(json.loads((ROOT / ".zenodo.json").read_text())["version"] == "1.1.0", ".zenodo.json version mismatch")

    # No bytecode or historical defective solver should enter the release tree.
    bad_bytecode = [p for p in ROOT.rglob("*") if p.is_file() and (p.suffix in {".pyc", ".pyo"} or "__pycache__" in p.parts)]
    require(not bad_bytecode, f"bytecode files present: {bad_bytecode[:5]}")
    bad_historical = [p for p in ROOT.rglob("*") if p.is_file() and ("HISTORICAL" in p.name or "historical_reference" in p.parts)]
    require(not bad_historical, f"superseded generalized-alpha files present: {bad_historical[:5]}")

    ga_root = ROOT / "reproducibility" / "generalized_alpha"
    single = rows(ga_root / "results" / "single_benchmark_corrected.csv")
    boundary = rows(ga_root / "results" / "boundary_groups" / "group_05.csv")
    compact_first = rows(ROOT / "data" / "section4_1_first_collision_corrected.csv")
    compact_second = rows(ROOT / "data" / "section4_1_second_collision_corrected.csv")
    require(single == compact_first, "first Section 4.1 collision table differs from production output")
    require(boundary == compact_second, "second Section 4.1 collision table differs from production output")
    require({r["newton_tol"] for r in single} == {"1e-10"}, "first collision tolerance mismatch")
    require({r["newton_tol"] for r in boundary} == {"1e-10"}, "second collision tolerance mismatch")
    require({r["newton_maxit"] for r in boundary if r["method"] == "GA"} == {"20"}, "second collision GA cap mismatch")
    require({r["newton_maxit"] for r in boundary if r["method"] == "M2"} == {"40"}, "second collision Pade cap mismatch")

    caps = rows(ROOT / "data" / "robustness_bank_newton_cap_audit.csv")
    total_ok = sum(int(r["status_ok"]) for r in caps)
    total_capped = sum(int(r["status_ok_with_any_cap_hit"]) for r in caps)
    require((total_capped, total_ok) == (2656, 4584), f"Newton-cap total mismatch: {total_capped}/{total_ok}")

    settings = {r["bank"]: r for r in rows(ROOT / "data" / "robustness_bank_settings_recovered.csv")}
    require(settings["bank01_mu_dt"]["dx"] == "0.05" and settings["bank01_mu_dt"]["c1"] == "1.25", "bank01 h/c1 mismatch")
    require(settings["bank02_dx_refinement"]["c1"] == "1.25", "bank02 c1 mismatch")
    require(settings["bank03_collision_severity"]["dx"] == "0.05", "bank03 h mismatch")
    require(settings["bank04_long_time"]["dx"] == "0.05" and settings["bank04_long_time"]["c1"] == "1.25", "bank04 h/c1 mismatch")

    required_figures = {
        "fig1_pg_symbol_filters_left.pdf",
        "fig1_pg_symbol_filters_right.pdf",
        "fig2_jacobian_h3_scaling_left.pdf",
        "fig2_jacobian_h3_scaling_right.pdf",
        "fig3_controlled_damping_tradeoff_left.pdf",
        "fig3_controlled_damping_tradeoff_right.pdf",
        "fig4_highthroughput_summary_left.pdf",
        "fig4_highthroughput_summary_right.pdf",
        "fig5_absolute_spatial_convergence_left.pdf",
        "fig5_absolute_spatial_convergence_right.pdf",
        "fig7_accuracy_cost_absolute_reference.pdf",
    }
    present = {p.name for p in (ROOT / "figures").glob("*.pdf")}
    require(required_figures <= present, f"missing current figure panels: {sorted(required_figures-present)}")

    print(f"SECTION4_COLLISION_TOLERANCE_OK first={len(single)} second={len(boundary)} tol=1e-10")
    print(f"ROBUSTNESS_NEWTON_CAP_TOTAL_OK {total_capped}/{total_ok}")
    print("RELEASE_METADATA_OK version=1.1.0")
    print("V1_1_0_RELEASE_AUDIT_OK")


if __name__ == "__main__":
    main()
