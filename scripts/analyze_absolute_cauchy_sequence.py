#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from analyze_absolute_fullfield_cauchy import load_npz, compare, write_csv


BASE = Path(__file__).resolve().parent

def find_core_pair(csv_path, n1, n2):
    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            if (
                r["kind"] == "grid_refinement"
                and int(r["ext_1"]) == 1
                and int(r["ext_2"]) == 1
                and int(r["N_1"]) == n1
                and int(r["N_2"]) == n2
            ):
                return float(r["E_core_eval12800"])
    return None


def apparent_order(e12, e23):
    if e12 <= 0 or e23 <= 0:
        return float("nan")
    return math.log(e12 / e23, 2.0)


def main():
    root = Path.cwd()
    p800 = root / "baseline_N800" / "cc_N800_ext1.npz"
    p1600 = root / "k22_ext" / "k22ext_cc_N1600_ext1.npz"
    p3200 = root / "k22_ext" / "k22ext_cc_N3200_ext1.npz"
    core_csv = root / "analysis" / "dweno_core_pairs.csv"

    missing = [p for p in (p800, p1600, p3200, core_csv) if not p.is_file()]
    if missing:
        print("MISSING_INPUTS")
        for p in missing:
            print(p)
        raise SystemExit(2)

    runs = {
        800: load_npz(p800),
        1600: load_npz(p1600),
        3200: load_npz(p3200),
    }

    out = root / "analysis_absolute_cauchy"
    out.mkdir(exist_ok=True)

    pair_results = {}
    for n1, n2 in ((800, 1600), (1600, 3200)):
        hist, selected = compare(runs[n1], runs[n2])
        prefix = out / f"absolute_{n1}_{n2}"
        write_csv(hist, str(prefix) + "_history.csv")
        core = find_core_pair(core_csv, n1, n2)
        last = hist[-1]
        pair_results[(n1, n2)] = dict(last=last, core=core, selected=selected)

        with open(str(prefix) + "_selected.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(selected[0].keys()))
            w.writeheader()
            w.writerows(selected)

    a = pair_results[(800, 1600)]
    b = pair_results[(1600, 3200)]

    metrics = [
        ("raw_relL2", "raw full field"),
        ("registered_relL2", "registered full field"),
        ("low_m16_registered", "registered low modes |m|<=16"),
        ("low_m32_registered", "registered low modes |m|<=32"),
        ("low_m64_registered", "registered low modes |m|<=64"),
        ("low_m128_registered", "registered low modes |m|<=128"),
        ("positive_registered_rel_to_fullfine", "positive-part difference/full-fine norm"),
        ("negative_registered_rel_to_fullfine", "negative-part difference/full-fine norm"),
    ]

    lines = []
    lines.append("ABSOLUTE DIRECT-WENO FULL-FIELD CAUCHY AUDIT")
    lines.append("Model: |K|(2,2), same-sign CC, ext=1, L=100, mu=0, T=10, CFL=.075, eps=1e-6")
    lines.append("")
    lines.append("FINAL-TIME PAIRS")
    for (n1, n2), d in pair_results.items():
        r = d["last"]
        lines.append(f"{n1}->{n2}")
        lines.append(f"  raw={r['raw_relL2']:.12e}")
        lines.append(f"  registered={r['registered_relL2']:.12e}  shift={r['global_shift']:+.12e}")
        lines.append(f"  low16={r['low_m16_registered']:.12e}")
        lines.append(f"  low32={r['low_m32_registered']:.12e}")
        lines.append(f"  low64={r['low_m64_registered']:.12e}")
        lines.append(f"  low128={r['low_m128_registered']:.12e}")
        lines.append(f"  positive/full={r['positive_registered_rel_to_fullfine']:.12e}")
        lines.append(f"  negative/full={r['negative_registered_rel_to_fullfine']:.12e}")
        lines.append(
            f"  negative_L2_fraction coarse={r['coarse_negative_L2_fraction']:.12e} "
            f"fine={r['fine_negative_L2_fraction']:.12e}"
        )
        lines.append(
            f"  min_u coarse={r['coarse_min_u']:.12e} fine={r['fine_min_u']:.12e}"
        )
        if d["core"] is not None:
            lines.append(f"  fitted-core pair difference={d['core']:.12e}")
            lines.append(
                f"  registered-full/core={r['registered_relL2']/d['core']:.6f}"
            )

    lines.append("")
    lines.append("CONTRACTION  (1600->3200)/(800->1600)")
    rows = []
    for key, label in metrics:
        e12 = a["last"][key]
        e23 = b["last"][key]
        ratio = e23 / e12
        p = apparent_order(e12, e23)
        lines.append(
            f"  {label}: E12={e12:.12e} E23={e23:.12e} "
            f"ratio={ratio:.9f} p_app={p:.6f}"
        )
        rows.append(
            dict(metric=key, label=label, E_800_1600=e12, E_1600_3200=e23,
                 contraction_ratio=ratio, apparent_order=p)
        )

    if a["core"] is not None and b["core"] is not None:
        ratio = b["core"] / a["core"]
        p = apparent_order(a["core"], b["core"])
        lines.append(
            f"  fitted core: E12={a['core']:.12e} E23={b['core']:.12e} "
            f"ratio={ratio:.9f} p_app={p:.6f}"
        )
        rows.append(
            dict(metric="fitted_core", label="fitted core",
                 E_800_1600=a["core"], E_1600_3200=b["core"],
                 contraction_ratio=ratio, apparent_order=p)
        )

    with open(out / "absolute_cauchy_contraction.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    lines.append("")
    lines.append("DECISION")
    reg_ratio = b["last"]["registered_relL2"] / a["last"]["registered_relL2"]
    low32_ratio = b["last"]["low_m32_registered"] / a["last"]["low_m32_registered"]
    if reg_ratio < 1 and low32_ratio < 1:
        lines.append(
            "  Both the registered full field and the fixed low-mode component contract "
            "from 800->1600 to 1600->3200."
        )
        lines.append(
            "  Therefore the previous statement that the absolute-flow full field is non-Cauchy "
            "at these resolutions is not supported by this refinement sequence."
        )
    else:
        lines.append(
            "  At least one of registered full-field or low-mode differences does not contract."
        )
        lines.append(
            "  The coherent cores are more reproducible than the full field; the section-5 "
            "core/remainder distinction remains necessary."
        )
    lines.append(
        "  This is a finite-grid refinement diagnostic, not a proof of continuum convergence."
    )

    report = "\n".join(lines) + "\n"
    (out / "absolute_cauchy_report.txt").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
