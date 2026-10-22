# v1.1.0

Reproducibility expansion and numerical-audit release.

- adds the corrected generalized-alpha implementation and validation;
- adds all 396 microbank rows, all 64 fine-boundary rows, the seven-run first
  Section 4.1 collision table, and the complete second-collision table;
- adds the Picasso scheduler scripts, manifests, logs and aggregate checks;
- adds collision-Jacobian and fixed-state mesh-scaling drivers, raw snapshots,
  frozen outputs and independent verification;
- documents that the historical Jacobian driver was not recovered byte-for-byte
  and that the included reconstruction reproduces all archived outputs;
- adds all four robustness-bank `cases.csv` files and complete per-run summaries;
- adds the MATLAB R2024b production chain, original case generators and SLURM
  launchers for the robustness screens;
- adds exact audits of observation windows, Newton settings, relative amplitude
  stopping criterion, and completed runs containing capped Newton solves;
- adds editable Python sources and CSV inputs for all manuscript figures;
- updates the one-command freeze check to cover all v1.1.0 components.

No superseded generalized-alpha result is used as valid data in this release.
The multi-gigabyte robustness-bank trajectory files are not duplicated; every
case can be rerun from the included production code and case table.

# v1.0.0

Initial public reproducibility release for the Direct-WENO and
Petrov-Galerkin coherent-core benchmark and invariant-balance checks.

- Default and full audits run in disposable copies and do not mutate the frozen tree.
