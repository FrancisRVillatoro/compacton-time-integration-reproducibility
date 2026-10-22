#!/usr/bin/env python3
"""One-command verification entry point for release v1.1.0."""
from pathlib import Path
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


EXPECTED_RAW = {
    "cc_N800_ext0.npz": "29babc260877594c9f041ffb05dd0853e0161facc0998d6e26a7057de4d7d921",
    "cc_N800_ext1.npz": "c214e367bba3e383b85c80881809f900899a34b7a30437662dc375d472be6a5d",
    "k22ext_cc_N1600_ext1.npz": "95603ba5531ae0fd90fdc5a9d08c1c54bba3f039b2db85691ed2cf27567b2eac",
    "k22ext_cc_N3200_ext0.npz": "050b4a9c64f9c10cd90df60765cfcbe457af3d40b59b7f85befe0ce2db302fc8",
    "k22ext_cc_N3200_ext1.npz": "a70fbfa816b97f07535fc203ca56a43a5e2b184618e51ee3e4f0ab08b0550560",
}


def run(path: Path, *args: str, cwd: Path | None = None) -> None:
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["MPLBACKEND"] = "Agg"
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    subprocess.check_call([sys.executable, str(path), *args], cwd=cwd or ROOT, env=env)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Verify the frozen v1.1.0 release. The default audit checks hashes, "
            "all v1.1.0 provenance chains and every figure. --full additionally "
            "recomputes the v1.0.0 Direct-WENO/Petrov-Galerkin benchmark from raw data."
        )
    )
    parser.add_argument(
        "--full", action="store_true",
        help="also rerun the expensive v1.0.0 raw benchmark chain",
    )
    parser.add_argument(
        "--_in-place", action="store_true", help=argparse.SUPPRESS
    )
    args = parser.parse_args()

    # Every audit regenerates some derived tables, reports, or PDFs.  Execute
    # it in a disposable copy so the frozen release tree remains byte-stable.
    if not args._in_place:
        with tempfile.TemporaryDirectory(prefix="compacton-v1.1.0-audit-") as tmp:
            sandbox = Path(tmp) / ROOT.name
            shutil.copytree(ROOT, sandbox)
            env = dict(os.environ)
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            env["MPLBACKEND"] = "Agg"
            for key in (
                "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
            ):
                env[key] = "1"
            command = [
                sys.executable,
                str(sandbox / "scripts" / "freeze_check.py"),
                "--_in-place",
            ]
            if args.full:
                command.append("--full")
            subprocess.check_call(command, cwd=sandbox, env=env)
        print(
            "FULL_AUDIT_DISPOSABLE_COPY_OK"
            if args.full else "FREEZE_AUDIT_DISPOSABLE_COPY_OK"
        )
        return

    for name, expected in EXPECTED_RAW.items():
        path = ROOT / "raw" / "direct_weno" / name
        got = sha(path)
        print(name, got, "OK" if got == expected else "FAIL")
        if got != expected:
            raise SystemExit(2)

    if args.full:
        for name in (
            "recompute_dweno_N3200_raw_checks.py",
            "recompute_absolute_spatial_from_raw.py",
            "recompute_table3_absolute_reference.py",
            "recompute_i2_table3.py",
            "build_final_table3.py",
            "make_fig5_spatial_convergence.py",
            "make_fig7_accuracy_cost.py",
        ):
            run(ROOT / "scripts" / name)
    else:
        print("V1_0_RAW_RECOMPUTATION_SKIPPED (use --full to enable)")

    # v1.1.0 additions and manuscript-facing artifacts.
    run(ROOT / "reproducibility" / "generalized_alpha" / "verify_generalized_alpha_traceability.py", "--full")
    run(ROOT / "reproducibility" / "jacobian" / "verify_jacobian_traceability.py")
    run(ROOT / "scripts" / "audit_robustness_banks.py")
    run(ROOT / "figure_sources" / "regenerate_and_verify.py", cwd=ROOT / "figure_sources")
    run(ROOT / "scripts" / "make_all_figures.py")

    # Nested scientific drivers may create bytecode while auditing.  Remove it
    # from the disposable audit copy before validating the release-tree hygiene.
    for cache_dir in ROOT.rglob("__pycache__"):
        shutil.rmtree(cache_dir, ignore_errors=True)
    for pattern in ("*.pyc", "*.pyo"):
        for bytecode in ROOT.rglob(pattern):
            bytecode.unlink(missing_ok=True)

    run(ROOT / "scripts" / "verify_v1_1_0_release.py")

    print("FREEZE_CHECK_V1_1_0_OK")


if __name__ == "__main__":
    main()
