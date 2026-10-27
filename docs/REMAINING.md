# Remaining repository work

Updated 26 October 2026. These are the outstanding code-integration and release
tasks; the supplied analysis workflows are described in the README and Coverage.

- [x] **P1 — IgFold generation code and settings.** Supplied in
  `kaggle/kernels/igfold-bench` with the retained coordinates, model-checkpoint
  hashes, software environment and per-target scores in
  `outputs/structure/igfold/generation`. `src/score_igfold_generation.py` rescores the
  retained coordinates and compares them with the primary record.
- [ ] **D1 — Docking evaluation procedure.** Add the DockQ and epitope-scoring
  script, notebook or exact commands. Identify the reference, chain mapping and
  scoring inputs, and provide the inputs or a usable retrieval/regeneration
  route. Prepared structures, regional scores and summary tables are supplied.
- [ ] **P3 — Licence and DOI release.** Agree the licence for original code,
  retain third-party terms, test the integrated workflows, refresh the package
  manifest and publish the declared version in a DOI-assigning repository.
  Add the version-specific DOI to the README, manuscript and response letter.

D1 completes the outstanding docking-evaluation component.
P3 completes the code-deposition requirement. The private GitHub collaboration
copy is not the DOI deposit.
