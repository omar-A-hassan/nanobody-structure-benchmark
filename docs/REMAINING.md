# Revision status and final inputs

Updated 26 October 2026. DONE means the author revision is complete by the
analysis, clarification or claim withdrawal stated below. It does not mean
that an unperformed experiment was completed or that the reviewers accepted it.

## Completed revisions

| Reviewer reference | Status | Completed work |
|---|---|---|
| R1 major 1; R2 comment 1, redundancy; R3 major 1a; R4 major 4, construction | DONE | Dataset selection, bound/free composition, identifiers, missing observations, sequence clustering and cohort membership are documented. |
| R1 major 2; R2 comment 1; R3 major 1b; R4 major 4, overlap | DONE: available-source analysis and scope response | Model-specific PDB/date/VHH/H3 screens and the 393-target sensitivity are supplied. The response explains why no antigen/epitope/interface screen is added. Neither VHH nor complex-level training independence is certified. |
| R1 minor, failures; R3 major 1c | DONE | Nine NanoBodyBuilder2 refusals and the separate partial output are identified and explained. |
| R2 comment 2; R3 major 2a/3c; R4 major 2 | DONE | All-loop benchmark metrics and all four prepared-input regional rows are in the paper. The prepared-input scorer, files, results and tests are supplied. |
| R3 major 2b | DONE: unsupported claim omitted | Verified SimpleFold/NanoBodyBuilder2 batch timings have defined boundaries. No IgFold runtime or three-model speed comparison is reported. No IgFold timing archive or new timing run is required. |
| R3 major 4a/5a/5b/6a; R4 minor 2 | DONE | Distributions, paired intervals, hard/easy definitions, length adjustment, probe labels and per-panel populations are supplied. |
| R1 major 3; R3 major 3a; R4 major 1/3/5 | DONE: claims narrowed | Docking is one illustrative holo re-docking case. H3-specific causation and general downstream-effect claims are withdrawn. A new contact study, repeated docking, second docking method or multi-complex benchmark is not promised. |
| R1 major 4 | DONE: scope-based response | The response explains the common VHH-only prediction task and why antigen-context prediction followed by docking is not added. It does not claim that comparison was performed. |
| R3 major 3b/7a/7b and minor 2 | DONE | Case selection, bounded architectural interpretation, joint-predictor literature and proportional conclusions are supplied. |
| R1 minor, Figure 1; R3 minor 1; R4 minor 1 | DONE | The workflow figure, docking-branch count and rounding are corrected. |
| R3 major 8a, structural/sampling analyses | DONE | Frozen inputs, scores, analysis code, sampling generators, cached probe data, run instructions and verification records are supplied. Full-paper code delivery still needs P1 and D1 below. |

## Only remaining handover inputs

- [ ] **P1 — IgFold generation code and settings.** Add the benchmark script
  or notebook, its model/weight version, refinement/numbering settings and
  execution environment. State a seed if one was used. Connect the supplied
  code to the frozen inputs and regional scorer. Fill manuscript P1 from
  those facts. No runtime measurement or full-benchmark rerun is requested.
- [ ] **D1 — Docking evaluation code and its data.** Add the script/notebook
  or commands that produce the reported DockQ and epitope metrics and table
  summaries. Include its reference/chain mapping, scoring inputs and per-pose
  values supporting the top-1 and best-of-top-10 results. Identify the supplied
  prepared files used by that workflow. There is no separate requirement for
  all 40 historical pose files, permanent server pages or a job-hash ledger.
- [ ] **P3 — Licence and deposited DOI.** Agree the licence for original code,
  retain third-party terms, publish the declared version in a DOI-assigning
  repository, and put its actual DOI in Code availability and the response.
  Check the affected workflows after P1/D1 are added and refresh the manifest.

P1/D1 complete the remaining full-paper code delivery under R3 major 8a.
P3 completes the editor's explicit deposit request. The private GitHub
repository is the collaboration copy, not the DOI deposit.

## Already supplied; do not request again

Sampling generation code, input sequences, sampler settings and saved feature/
result tables are present; see [Sampling](SAMPLING.md). No complete trajectory
archive, compulsory sampling GPU rerun, historical timing CSV, redistributed
model weights or archive of unreported experiments is required. Preserve
existing records. Do not describe saved-data replay as a new GPU experiment.
