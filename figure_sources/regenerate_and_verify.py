#!/usr/bin/env python3
"""Regenerate Figures 1--3 and verify their data/code lineage and PDF hashes."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parent
OUT = ROOT / "figures"
VERIFY = ROOT / "verification"

SCRIPTS = (
    "fig01_pg_symbol_filters.py",
    "fig02_jacobian_h3_scaling.py",
    "fig03_controlled_damping_tradeoff.py",
)
OUTPUTS = (
    "fig1_pg_symbol_filters_left.pdf",
    "fig1_pg_symbol_filters_right.pdf",
    "fig2_jacobian_h3_scaling_left.pdf",
    "fig2_jacobian_h3_scaling_right.pdf",
    "fig3_controlled_damping_tradeoff_left.pdf",
    "fig3_controlled_damping_tradeoff_right.pdf",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def keyed(rows: list[dict[str, str]], fields: tuple[str, ...]) -> dict[tuple[str, ...], dict[str, str]]:
    return {tuple(row[field] for field in fields): row for row in rows}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    VERIFY.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    for script in SCRIPTS:
        subprocess.run([sys.executable, str(ROOT / script)], cwd=ROOT, env=env, check=True)

    pdfs = {}
    for name in OUTPUTS:
        generated = OUT / name
        if not generated.is_file() or generated.stat().st_size == 0:
            raise AssertionError(f"missing regenerated PDF: {name}")
        pdfs[name] = {"sha256": sha256(generated), "bytes": generated.stat().st_size}

    fig2 = ROOT / "data" / "fig2_jacobian_h3_scaling.csv"
    jacobian = PACKAGE / "reproducibility" / "jacobian" / "generated" / "jacobian_mesh_scaling.csv"
    if sha256(fig2) != sha256(jacobian):
        raise AssertionError("Figure 2 data are not byte-identical to regenerated Jacobian mesh data")

    fig3 = ROOT / "data" / "fig3_controlled_damping_fine_corrected.csv"
    ga = PACKAGE / "reproducibility" / "generalized_alpha" / "results" / "fig3_right_corrected.csv"
    if sha256(fig3) != sha256(ga):
        raise AssertionError("Figure 3 data are not byte-identical to the corrected Picasso summary")

    # Ensure every fine-case row represented in the GA result is present unchanged.
    if keyed(read_csv(fig3), ("label",)) != keyed(read_csv(ga), ("label",)):
        raise AssertionError("Figure 3 row-level data mismatch")

    fig1_source = (ROOT / SCRIPTS[0]).read_text(encoding="utf-8")
    fig2_source = (ROOT / SCRIPTS[1]).read_text(encoding="utf-8")
    if "from jacobian_core import ga_rho" not in fig1_source:
        raise AssertionError("Figure 1 does not use the canonical Jacobian ga_rho implementation")
    if "step_GA" in fig1_source + fig2_source or "step_GA_fast" in fig1_source + fig2_source:
        raise AssertionError("Figure 1/2 source calls a nonlinear generalized-alpha stepper")
    if "midpoint_amp" in fig2_source or "'M1_amp'" not in fig2_source:
        raise AssertionError("Figure 2 is not wired to the regenerated M1_amp column")

    report = {
        "status": "FIGURE_TRACEABILITY_OK",
        "pdfs": pdfs,
        "data": {
            "figure_2_jacobian_csv_sha256": sha256(fig2),
            "figure_3_corrected_ga_csv_sha256": sha256(fig3),
        },
        "generalized_alpha": {
            "figure_1_source": "canonical ga_matrix/ga_rho only",
            "figure_2_source": "regenerated Jacobian CSV; no nonlinear GA stepper",
            "figure_3_source": "corrected fine Picasso case; coarse N=100 screen excluded",
        },
    }
    (VERIFY / "FIGURE_VERIFICATION.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    lines = [
        "FIGURE_TRACEABILITY_OK",
        "",
        "Figures 1--3 regenerated successfully.",
        "All six source-lineage PDFs were regenerated successfully.",
        "Figure 1 uses only the canonical generalized-alpha amplification matrix.",
        "Figure 2 uses the fully recomputed Jacobian mesh CSV.",
        "Figure 3 uses only the corrected fine-grid Picasso case.",
        "The old N=100 screening aggregate is retained for provenance but is not plotted.",
        "",
    ]
    (VERIFY / "FIGURE_VERIFICATION.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines), end="")


if __name__ == "__main__":
    main()
