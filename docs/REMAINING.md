# Remaining work for manuscript reproducibility

The saved structural analyses and the SimpleFold/NanoBodyBuilder2 inference
repeats are included. Resolve the following items before describing this
package as the complete code and data archive for the manuscript. Check an
item only after its files and verification are present.

## IgFold inference and runtime

- [ ] Add the code or notebook that generated the 597 benchmark predictions.
- [ ] Record the executed version, weights or hashes, seed, refinement and
  numbering settings, environment and hardware.
- [ ] Recover `prediction_results.csv` or the original timing log. Verify the
  reported 4.08 seconds per target against its timing boundary, output count
  and parallel execution before making a speed comparison.
- [ ] Add the prediction coordinates and replay the regional scorer. Version
  any new predictions separately and compare them with the primary results.

The package contains all 597 IgFold regional scores, six visualization
examples and an IgFold code path in the separate TNP notebook. That TNP run
built zero IgFold structures; it is not the benchmark generator. The evidence
is recorded in `provenance/igfold_benchmark_evidence.json`.

Addresses R3 major 2b and 8a. Manuscript author input P1 marks the article
details that require these records.

## Docking case

- [ ] Add the four exact prepared VHH inputs, fixed receptor, experimental
  complex, file hashes and residue/chain correspondence.
- [ ] Record the case-selection reason and preparation code/settings.
- [ ] Score H1, H2, H3 and framework for the prepared inputs, including the
  prepared experimental control. Do not substitute unverified benchmark rows.
- [ ] Add HDOCK settings and job records, the retained top-ten poses for each
  input, and the pose-ranking rule.
- [ ] Add DockQ version/commands, chain mapping, native-contact and epitope
  definitions, per-pose outputs, and table-generation code.
- [ ] Reproduce every docking and epitope value retained in the manuscript.

Addresses the 7R73 part of R2 comment 2, R3 major 3b/3c/8a and R4 major 2.
Manuscript author inputs P2a and P2b mark the required article content.

## Other reviewer requests

- [ ] Resolve R1 major 4 with the requested antigen-context comparison or an
  author-approved, justified scope response. No such experiment is included.
- [ ] Resolve the antigen/epitope/interface component of R1 major 2 separately
  from the VHH overlap screen.
- [ ] Resolve R4 major 3 and 5 with regional native-interface evidence and
  docking controls/repeats, or an explicit justified response supported by the
  final case-study scope. Scope changes are not experiments.

The available-source training audit has incomplete manifest and pretraining
coverage. The 393-target subset is not certified training-independent. Any
extended audit must preserve its sources, thresholds and new cohort lists.

## Sampling inputs

- [ ] Complete raw sample/trajectory coordinate archives and their run records
  where required for end-to-end regeneration. Saved summaries and supplied
  probe embeddings support the documented analysis replay.

## Public release

- [ ] Confirm redistribution rights and select a licence for original code.
  Preserve third-party data and code terms.
- [ ] Add `LICENSE` and `CITATION.cff` with verified creator and version data.
- [ ] Run analysis replay, coordinate rescoring and tests on the final release.
- [ ] Publish the versioned archive in a DOI-assigning repository and verify
  its public files. A draft or reserved DOI is not a published deposit.
- [ ] Insert the version DOI in Code availability and the response letter.
  Confirm that declared coverage matches the deposited files.

Addresses the editor's code-deposit requirement and R3 major 8a. Manuscript
author input P3 marks the DOI and release-coverage check.
