# Entry points and dependency graph

The table lists the analysis entry points, their inputs, outputs and execution
environments. Structural analyses and reference reconstruction use repository
data and recorded upstream sources; manuscript assembly is outside their scope.

| Entry point | Inputs | Output | Environment |
|---|---|---|---|
| `reviewer_revision_analysis.py` | Frozen targets; three original score JSONs; pairwise identity cache or sequences; hard, trajectory and hit-rate results; crystal labels; NBB manifest | Regional statistics/CSVs, redundancy and uncertainty statistics, distributions | Analysis |
| `revision_training_audit.py` | Complete compressed SAbDab2 table; NBB manifest; derived IgFold VH/H3 pool; nearest-training caches; pairwise identities; original scores and release dates | Initial 398-representative lower-identified-risk analysis and overlap tables | Analysis |
| `revision_overlap_extension.py self-test` | Saved source-bound sequence pools | Pruned-vs-brute-force alignment test and H3-span regression | Overlap |
| `revision_overlap_extension.py summarize` | Initial audit, three expanded nearest-neighbour tables, pairwise matrix, original scores | Expanded flags, 393 representatives, cohort scores and paired/PDB-block sensitivity | Overlap |
| `revision_overlap_extension.py analyze-simplefold` | Public SwissProt/AFESM IDs, extended dictionary, repository/recipe evidence, benchmark entity mapping | Identifier/recipe scope diagnostics | Overlap |
| Other overlap collection stages | Public sources and receipts; ANARCII for domain numbering | Entity, sequence-pool and nearest-neighbour evidence | Collection |
| `reconstruct_revision_references.py` | Source release, SD-H table, frozen targets, scores and cohort CSVs | Exact frozen replay, exclusions, alternate-atom diagnosis and policy displacement | Reference |
| `revision_reference_sensitivity.py` | Validated reconstruction summary, 393-cohort CSV, initial audit and original scores | Original-score subsets excluding affected references | Analysis |
| `kaggle/kernels/sf-traj-tfm-597/` | Frozen 597 inputs, SimpleFold/ESM weights and the notebook's sampler logging patch | Six generated samples per target, intermediate per-residue coordinate averages and `traj_results.json` | Supplied GPU source; settings and instructions in `docs/SAMPLING.md`; separate from saved-data replay |
| `kaggle/kernels/sf-hard/` | Frozen inputs and the explicit 36-target hard/easy lists | Six-sample generation, aligned spread/errors and `hard_results.json` | Historical GPU notebook; separate from saved-data replay |
| `kaggle/kernels/sf-hitrate/` | Frozen inputs, 30 hard/ten easy lists and four sampling settings | 3,200 generated predictions and per-setting `hitrate_results.json` | Historical GPU notebook; separate from saved-data replay |
| `src/tnp/tnp_crystals.py` | Frozen references and `tnp_metrics.py` | Crystal conformation labels; percentile-calibrated geometry | Analysis |
| `esm_probe_analyze.py` | Full retained embeddings, 518 labels/IDs, framework sequences | All-layer linear/MLP probe and shuffled controls | Analysis |
| `reproduce_esm_probe.py` | Same embeddings/labels and unchanged original fitting source | Results, 110 sets of out-of-fold probabilities/labels/groups and environment receipt | Analysis |
| `recreate_reused_figures.py` | Unchanged original figure source, trajectory/hit-rate/ESM saved data | Three reusable figure pairs; average-tie Spearman validation | Analysis |
| `make_revision_workflow.py` | Verified common/expanded cohort counts | Study workflow diagram; docking branch is descriptive only | Analysis |
| `verify_retained_structure_scoring.py` | Legacy score functions, frozen references, separate visualization PDBs/reports | Limited retained-example scoring checks | Analysis |
| `score_prepared_vhh.py` | Four prepared VHH PDBs, frozen 7R73A reference and original Kabsch source | Regional RMSD CSV/JSON with input hashes and independent alignment checks | Analysis; `--check` verifies saved results without replacing them |
| `audit_complex_overlap.py` | Source antigen annotations, model-source PDB identifiers and supplied exact-case source coordinates | Antigen sequence matches, native heavy-domain contacts, epitope/contact-pattern comparisons and per-target summaries | Analysis; see `docs/COMPLEX_OVERLAP.md` |
| `kaggle/kernels/7r73-docking/docking.ipynb` | Frozen VHH sequence, downloaded 7R73, supplied IgFold PDB, SimpleFold/NBB2 sources and weights | Prepared VHH and antigen PDBs | Supplied execution; moving source/dependencies, no new regeneration test; no downstream pose scoring |
| `replay_revision_nbb2.py` | Three selected VHHs, NBB2/numbering/OpenMM stack and separately acquired weights | Optional new CPU predictions and provenance, distinct from primary benchmark | Existing historical replay lock; not run during assembly |
| `verify_archive_reruns.py` | Versioned completed repeat ZIP/member manifests and scores, frozen targets, recovered original scorer and primary scores | Archive integrity, exact v1 coordinate/scorer replay, checkpoint logging status and original/repeat differences | Optional inference-check; NBB2 v1 locally verified |
| `prepare_archive_reruns.py` | Recovered original inference notebooks; branch revision resolved when preparing SF | Coordinate-retaining notebook source and metadata; corrected extensionless NBB2 weight logging | Executed notebooks and completed coordinate archives have separate version receipts |
| Other `src/tnp/` tools | Original saved result collections | Historical conformation summaries | Analysis; revised statistics take precedence |

The original notebooks retain model installation/inference instructions and
their original execution assumptions. Analysis locks do not promise to resolve
the original GPU inference software. Source hashes and historical execution
limits are recorded under `provenance/`.

Sampling inputs, settings, execution order and outputs are documented in
`docs/SAMPLING.md`. Intermediate-frame archives are optional. Saved-result
analysis replay does not claim to rerun GPU generation.

`provenance/original_scripts/consolidate_stats.py` is an unchanged historical
aggregator, outside the active analysis entry points. It also expects unrelated
tabular diagnostic inputs that are not part of this structural package; it is
retained for provenance only. The manuscript builders remain excluded.

The versioned repeat verifier uses temporary coordinate extraction and an
explicit separate report path. It preserves the imported v1 verification,
accepts the explicitly missing v1 checkpoint map, and can require four complete
`nanobody_model_1`–`nanobody_model_4` hash entries for newer runs. Logged hash
values are provenance; model weights are neither shipped nor rehashed locally.

## Public package adaptations

- The complete 15,641-row `ab_split.csv` is losslessly gzipped, preserving its
  original decompressed checksum. The public audit reads gzip; reconstruction
  validates the corresponding extracted archive member against its gzip bytes.
- The IgFold experimental training ZIP is replaced by its exact derived 2,649
  unique heavy/H3 pairs and 4,275 PDB identifiers. The original extraction block
  is retained in `provenance/original_scripts/revision_training_audit.py` and
  the runnable `tools/extract_igfold_sequence_pool.py`. The source archive MD5,
  SHA256 and public retrieval URL are included. Original source hash fields in
  numerical audit results continue to identify the original source bytes.
  One unmodified 7NXB archive member is retained under `data/regression/` so
  the extension self-test exercises the original H3-span regression offline.
- The large SimpleFold extended identifier dictionary is losslessly gzipped.
  The corresponding analysis reads the compressed copy. It remains an
  identifier dictionary for the extended recipe, not a complete training pool.
- The full ESM arrays are sliced along layers into five NPZ files, with dtype,
  shape and C-order array checksums. Reassembly checks every array value; the
  original compressed NPZ is an ignored local reconstruction, avoiding a
  272 MB public Git blob.
- Reused-figure selection resolves its source within this repository. The
  original `make_figures.py` is unchanged, including its historical tie helper
  and superseded benchmark bar plot. The portable driver calls only the three
  reusable figures and supplies the revised average-tie Spearman definition.
- Probe fitting and result fields are unchanged. The public probe source uses
  descriptive comments and exploratory diagnostic labels; the original selected
  source is retained in `provenance/original_scripts/esm_probe_analyze.py`.
  Notebook Markdown descriptions are clarified without changing code cells.
- The public overlap receipt explicitly covers the copied subset, including
  compressed sources and the adapted audit. It verifies those entries before
  and after each extension command. Replay should start from a fresh copy;
  rerunning earlier analyses with changed dependencies can invalidate it.

Numerical transformations, cohort definitions, bootstrap seeds, thresholds,
original saved scores and frozen reference coordinates were not altered.
Hash-bearing transitive locks are under `environments/`. `uv.lock` locks the
Python 3.13 project; Python 3.12 commands intentionally use `--no-project`.
