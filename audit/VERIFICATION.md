# Verification summary for v1.1.0

Run from the repository root:

```bash
python scripts/freeze_check.py
```

A successful verification ends with:

```text
FREEZE_CHECK_V1_1_0_OK
```

The default check validates the frozen v1.0.0 raw-input hashes and the complete
v1.1.0 provenance chain. Running `python scripts/freeze_check.py --full` also
recomputes the v1.0.0 Direct-WENO/Petrov-Galerkin benchmark from raw trajectories.
The v1.1.0 audit:

1. reruns the corrected generalized-alpha matrix, midpoint and order tests;
2. reconstructs the generalized-alpha aggregates from the 44 Picasso group CSVs;
3. checks all retained scheduler logs and the two Section 4.1 collision tables;
4. compares regenerated Jacobian spectra and mesh-scaling fits with the frozen
   numerical outputs;
5. audits every robustness-bank case and summary row, the MATLAB stopping rules,
   the relative factor-50 amplitude criterion and the Newton-cap counts;
6. regenerates every manuscript figure panel from its frozen CSV source.

The full collision trajectory and ARPACK Jacobian experiments can be rerun with:

```bash
python reproducibility/jacobian/verify_jacobian_traceability.py --full
```

The default freeze check and both full commands are non-mutating: all derived
outputs are created and checked in disposable copies, leaving the frozen release tree unchanged.

The historical Jacobian eigensolver driver was not recovered byte-for-byte;
the reconstructed driver and this qualification are part of the frozen release.
