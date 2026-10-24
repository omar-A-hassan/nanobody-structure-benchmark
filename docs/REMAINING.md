# Remaining work for manuscript reproducibility

Updated 24 October 2026. This checklist distinguishes workflow delivery,
claim-dependent evidence, scientific requests and optional archival work.
An unchecked item is not automatically a request for a new full GPU run.
Check an item only after the stated evidence is present. This document does
not certify reviewer acceptance or change the manuscript.

## 1. IgFold workflow

Already included: all 597 regional-score records, six visualization examples
and an IgFold code path in a separate TNP notebook. That TNP run built zero
IgFold structures; it is not the benchmark generator. See
`provenance/igfold_benchmark_evidence.json`.

- [ ] Supply the original benchmark generation procedure, or a documented
  and validated regeneration procedure for the same 597 inputs. A script or
  notebook is sufficient; no particular historical filename is required.
- [ ] Record the code version, model/weight identifier and retrieval source,
  refinement and numbering settings, seed where used, and tested environment.
  Do not invent unknown historical versions or seeds.
- [ ] Connect generation to the supplied regional scorer. Use original
  prediction files if recovered. Otherwise retain regenerated outputs as a
  separate version, compare per-target scores with the primary records, and
  report differences. New predictions are not recovered original predictions.

**Runtime is claim-dependent.** If the paper retains an IgFold runtime or
speed comparison, provide measurement records, timing boundaries, output
count, hardware and concurrency. `prediction_results.csv` is one possible
source, not a required filename. The historical 4.08 s/target report is not
verified timing evidence. If adequate records cannot be obtained, the authors
must remove that claim or report a separately measured, labelled repeat.
Do not rerun the entire benchmark solely to recreate the original CSV name.

Addresses R3 major 2b/8a and manuscript input P1. Available model identifiers
and retrieval instructions can identify weights without redistributing them.

## 2. Docking workflow and reported inputs

- [ ] Supply the exact four prepared VHH files, fixed antigen, experimental
  complex and residue/chain correspondence used for the reported case. Record
  file identities, preferably with checksums. Unverified benchmark predictions
  are not substitutes for these prepared inputs.
- [ ] Supply the preparation script or sufficiently detailed preparation
  steps, software versions and settings. Record HDOCK mode, restraints,
  receptor/ligand assignment and pose-ranking rule.
- [ ] Retain the poses needed to check the published endpoints. The current
  paper reports top-1 and best-of-top-10 results, so supply the retained ten
  candidates per input and their scores. Every unreported pose explored by
  the server is not required.
- [ ] Supply runnable regional, DockQ and epitope scoring and table-generation
  code, with versions, chain mapping, atom/contact definitions and per-pose
  results. Check the reported table values from the supplied prepared files
  and poses, including H1/H2/H3/framework scores for the prepared control.

A web-server run can be documented as manual submission followed by local
analysis. A new HDOCK API client or replacement server is not required. Keep
job IDs and original server reports if available; permanent access to an old
job page is not a substitute for, or a prerequisite to, retaining the inputs,
settings and outputs needed for the reported analysis.

Addresses the 7R73 part of R2 comment 2, R3 major 3c/8a and R4 major 2.
Inputs P2a/P2b also require the case-selection explanation and prepared-input
metrics in the manuscript. Case selection is an author response, not a new
software component.

## 3. Sampling generation and description

Sampling code is already included:

- [x] [Full-benchmark trajectory notebook](../kaggle/kernels/sf-traj-tfm-597/sf-traj-tfm-597.ipynb):
  six predictions for each of 597 inputs, 500 generation steps, logging every
  five steps, framework alignment and feature calculation.
- [x] [Hard/easy notebook](../kaggle/kernels/sf-hard/sf-hard.ipynb):
  target lists, six-sample generation and uncertainty scoring on 36 targets.
- [x] [Hit-rate notebook](../kaggle/kernels/sf-hitrate/sf-hitrate.ipynb):
  30 difficult targets and ten controls, 20 samples at each of four settings.
- [x] Saved per-target results and cached probe embeddings support the
  supplied analysis replay. The existing 15-job replay and five tests passed;
  those checks did not regenerate the GPU sampling experiments.

Samples are repeated SimpleFold predictions. Trajectories are intermediate
generation coordinates, not molecular-dynamics or docking trajectories.

- [ ] Specify fixed, usable source/dependency versions and model retrieval
  instructions for the included generators. Their current installation cells
  use a moving SimpleFold branch and some unpinned dependencies. A newly
  selected version must not be labelled as a recovered historical version.
- [ ] Document the run order, input attachment, GPU setup, output locations
  and expected target/sample counts. Verify the documented path before calling
  regeneration tested. Start with a representative GPU smoke check; a smoke
  check validates integration, not replication of full-cohort statistics.
  Full sampling repeats are needed only for a scientific validation question
  or a changed result, not merely to populate an archive.
- [ ] Resolve the trajectory-coordinate description. The included logger
  averages the unmasked atom coordinates within each residue; the manuscript
  Methods and Figure 3 currently call the logged values C-alpha coordinates.
  Confirm the executed logger and align the description with the actual
  feature. Do not silently relabel or replace numerical results.

**Not required by default:** every generated structure or every intermediate
frame. Full raw archives would allow an independent audit of the exact
historical coordinate-to-feature calculations. Their absence limits that
audit, but does not mean the generation code is absent or invalidate saved-data
analysis replay. Preserve available raw files; do not rerun all sampling
solely to recreate them. No reviewer explicitly requested a full frame archive.

## 4. Final code deposit

- [ ] Confirm redistribution rights and approve the licence/terms for original
  code. Supply those terms, normally in `LICENSE`, and preserve third-party
  notices. Upstream data and weights may remain at stable retrieval sources
  with identifiers and instructions; copying every upstream archive is not
  a blanket requirement.
- [ ] After integrating the missing workflows, test the affected entry points,
  scoring and tables from a clean checkout. Run the existing checks on the
  final package and update its manifest. Do not rerun every historical GPU
  experiment for a documentation-only change. Check scientific agreement;
  bit-for-bit equality across different GPU environments is not a blanket
  requirement. Preserve differences rather than silently replacing results.
- [ ] Publish the declared code/data contribution in a versioned DOI-assigning
  repository with verified creator, version and source metadata. Verify the
  accessible deposited files and insert the DOI in Code availability and the
  response. The editor explicitly requires this deposit. A private GitHub
  release or a reserved DOI does not complete it.

`CITATION.cff` is recommended for citation convenience, not a separate reviewer
requirement or a prerequisite to testing the code. The citation metadata in
the deposit must still be accurate. See [GitHub's citation-file documentation](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files).
This section addresses the editor's requirement, R3 major 8a and input P3.

## 5. Scientific responses, separate from packaging

The authors must resolve the following requests in the manuscript and response.
They are not automatically mandatory repository additions or instructions to
launch experiments. If an experiment is undertaken and reported, include the
inputs, procedure and outputs needed to check it. If it is not undertaken,
give an evidence-based justification and the corresponding scope/claim changes;
a justification does not guarantee reviewer acceptance.

- [ ] R1 major 4: antigen-context prediction and downstream comparison, or a
  justified response to that request. No such experiment is included.
- [ ] R1 major 2: antigen/epitope/interface overlap, distinct from the completed
  available-source VHH screen.
- [ ] R4 major 3: evidence about CDR-H3 contributions to the native interface.
- [ ] R4 major 5: docking controls/repeats and variance, or a justified,
  explicitly limited response. Case-study wording is not a new control run.
- [ ] R3 major 3b: explain case selection in the manuscript and response.

Full training independence remains unestablished. The 393-target set must not
be called leakage-free. This is not an instruction to recover every unavailable
pretraining corpus; it requires accurate coverage, limits and claims. Additional
screens, if performed, must retain their sources, thresholds and cohort lists.
