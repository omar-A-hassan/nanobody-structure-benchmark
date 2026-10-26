# Reproduction coverage

[Remaining work](REMAINING.md) lists outstanding integration and release tasks.
`reproduce.py` is the public entry point for
integrity checks and saved-analysis replay.

## Included scientific material

The 597 frozen targets retain all metadata, VHH sequences, three IMGT CDR
spans and the coordinate/null representation. The three saved score collections
contain 597 SimpleFold, 597 IgFold and 588 NanoBodyBuilder2 records, giving 588
common targets with eight regional/alignment scores each. Full and framework
fits, paired target bootstrap, PDB-block intervals, sequence redundancy,
available-source training overlap, reference policy and length-adjusted
uncertainty analyses are included.

The complex-source audit covers polymer-antigen sequences for 281 common-set
references and 280 usable native interfaces, with model-linked antigen,
epitope and heavy-domain contact-pattern matches. It records unavailable
partners and zero-occupancy exclusions separately. These measurements and
their limits are in `docs/COMPLEX_OVERLAP.md` and Supplementary Table S3.

The initial overlap audit has 398 nonredundant lower-identified-risk
representatives. The expanded available-source screen yields 393 from 362 PDB
entries. NBB2 recovery covers all 247 previously unmapped manifest PDBs, with
290 accepted heavy-domain candidates and 213 unique VH/H3 pairs. Candidate
recovery does not prove which chain trained a historical checkpoint.
The IgFold OAS pools screen 16,141 paired and 20,052 unpaired heavy records;
light-only FASTAs and resource-fork members are excluded. All heavy records
have an H3 assignment. Invalid coordinate members are explicitly diagnosed,
with valid FASTA/H3 information retained in the conservative sequence screen.

SimpleFold's retained public files contain 269,003 SwissProt and 2,036,414
AFESM identifiers. The extended dictionary has 8,616,595 member occurrences
for the extended recipe. Identifier matches and current entity/accession
checks do not establish complete VHH homology or historical checkpoint
membership. This repository does not certify training independence.

The standalone reference reconstruction selects 604 representatives from
3,414 single-domain rows, 3,350 SD-H rows and 2,404 eligible post-cutoff rows.
It reproduces all 597 frozen records, including 1,791 CDR spans, 72,378 residues,
843,498 coordinate components and 8,346 null atom slots. The seven other
representatives fail the explicit reconstructed checks. These checks are a
validated reconstruction, not recovered historical extraction code.

The frozen reference uses the last source row for duplicate alternate-location
atoms. Highest-occupancy selection changes 111/597 references, including
108/588 common targets. The supplied sensitivity subsets retain original
prediction scores; they do not rescore predictions against altered references.
The largest reference-to-reference framework-aligned H3 displacement is
approximately 0.6415 Å.

Saved sampling material includes the full 597-target trajectory summary,
available partial-run summaries/trajectory arrays, hard/easy uncertainty
inputs, and finite coverage for 30 difficult plus ten control targets across
four settings with 20 samples each. The retained artifacts do not constitute
all raw coordinates for every trajectory or sample. The original sampling,
trajectory, hit-rate, conformation and visualization notebooks are included.
Generation code, inputs, settings and saved analysis records are supplied;
[Sampling](SAMPLING.md) gives their use. Local replay evaluates saved-data
analyses. GPU generation and the notebooks' moving-source installations are
separate workflows.
The trajectory logger records per-residue means over unmasked atom coordinates,
not C-alpha-only coordinates. The Methods and Figure 3 caption use that definition.

The ESM probe includes all embeddings for 518 labelled targets (330 kinked,
188 extended; the other 79 reference labels are excluded), all 37 layers,
mean/anchor readouts, grouped five-fold fits, linear and MLP models and 50/20
shuffled controls. The complete rerun captures 110 cross-validation calls
(550 fitted folds). Best-layer selection is exploratory and is not nested
cross-validation. Threshold flags in the saved results are exploratory
diagnostics, not evidence of preregistration.
The exact historical scikit-learn version was not recorded; the validated
current rerun uses 1.9.1. Cached-embedding fitting is distinct from generating
new ESM2 embeddings with model weights.
The probe evaluates coarse-class decodability. It does not establish the
mechanism of SimpleFold coordinate errors. The original probe source is retained
under `provenance/original_scripts/`; fitting and stored results are unchanged.

## Inference provenance and unavailable work

- The docking preparation notebook and five prepared PDBs are included.
  All four VHH sequences match the 126-residue frozen 7R73A reference.
  Regional rescoring uses 492 matched N/CA/C/CB atoms and a 350-atom framework
  fit for each input. The CSV/JSON includes file hashes and an independent
  BioPython alignment check. The prepared control's computed regional RMSDs
  round to zero. These are monomer scores, not docked-pose scores.
  The two supplied docking CSVs are summary records, not a local pose-scoring
  replay. The remaining docking input is the evaluation procedure and its
  scoring inputs; see `docs/REMAINING.md`.

- The IgFold benchmark generation procedure and settings are not included.
  All 597 regional scores are included. The docking preparation notebook copies
  one existing IgFold PDB rather than generating the benchmark predictions.
  The paper reports no IgFold runtime or three-model speed comparison.
  The historical evidence record remains in `provenance/` for reference.
- Original SimpleFold/NBB2 full-benchmark notebooks explicitly deleted
  predicted coordinates after scoring. Original results and recovered source
  notebooks survive. Historical SimpleFold commit/checkpoint hashes and full
  filesystem iteration order, and the exact historical NBB2 installation,
  cannot be reconstructed from those saved artifacts.
- Visualization PDBs are included and checked as a limited scoring exercise.
  The historical IgFold record describes its six examples as derived from
  existing benchmark predictions; the displayed H3 scores agree with the
  primary records. These examples are not the full prediction batch or
  verified prepared docking inputs. SimpleFold/NBB2 visualization runs remain
  separate from their primary benchmark records.
- An existing three-target NBB2 CPU check retains three new predictions,
  scores, source/weight hashes and execution receipts. Its maximum absolute
  difference across 24 regional scores is 0.058 Å. It is a later pipeline check,
  not recovery of the original predictions.
- NBB2 version 1 is a separate completed repeat: 588 retained
  coordinates in one ZIP, manifest, original run verification, scores/logs,
  environment, run provenance, notebook and separate comparison JSON/CSV.
  All 588 primary IDs recur; local original-scorer replay matches the complete
  saved repeat rows exactly. Whole-chain H3 median RMSD is 2.4835 Å primary
  versus 2.474 Å repeat; median absolute change is 0.018 Å, maximum 1.967 Å.
  The nine unscored inputs are the same nine absent from the primary NBB2
  score set. This repeat neither replaces primary scores nor establishes
  identity with deleted historical coordinates.
- V1's checkpoint hash map remains empty. Its logging filter omitted
  extensionless `nanobody_model_*` weights; no weight hashes are inferred from
  other runs. The current preparation source includes corrected logging and
  a four-weight assertion. The completed NBB2 v2 and SimpleFold v1 runs are
  included as separate versions. Their local scorer replays match all
  saved rows exactly. NBB2 v2 records four weight hashes; SimpleFold records
  its 100M flow-checkpoint hash, source commit and worker input order. The
  ESM2 model identifier is recorded, but its weight-file hash was not logged.
  SimpleFold regional scores differ from the primary records by at most
  0.001 Å across 597 targets. NBB2 v2 retains the same 588 target IDs and nine
  refusals; its whole-chain H3 median is 2.474 Å, median absolute per-target
  change 0.014 Å and maximum change 1.952 Å.
  [Run-specific facts and local commands](NBB2_VERSION1.md) distinguish the
  archived v1 notebook from the newer preparation source.

## Deliberate exclusions and verification limits

Local working documents, manuscript builders, credentials, pretrained weights,
native executables, complete multi-GB source archives and incomplete downloads
are excluded. The source selection and
storage transformations are explicit in `provenance/source_selection.json`.
Full source archives remain available from their recorded upstream locations.
Their complete remote downloads and full ANARCII recollection are not part of
the saved-data verification unless explicitly run in the separate environment.
The auxiliary crystal-label replay gives identical labels, rounded metrics
and summary values. Unrounded alpha/tau angles differ by at most
2.842170943040401e-14 degrees, within the stated 1e-12-degree absolute check.
This is recorded as roundoff agreement rather than exact JSON equality.
The completed NBB2 v1, NBB2 v2 and SimpleFold v1 replays have separate reports.
Checkpoint/archive boundary tests reject corrupt hashes, duplicate members and
unsafe paths. The static inventory validates all completed repeat ZIPs and
their manifests. `verification/reexecution.json` records the analysis replay;
`verification/source_integrity.json` records comparison with the selected
sources. Earlier snapshot checksums are retained under `provenance/core_snapshot/`.

The verifier checks syntax, configured credential patterns, notebook output
removal, repository blob sizes, frozen targets/scores, curated receipt hashes
and, when sealed, the complete package manifest. Credential-pattern scanning
has finite coverage. These verifier checks do not publish a repository, DOI
record or Kaggle dataset. A separate private GitHub collaboration release
exists; it is not a public DOI deposit. No project licence is inferred from
third-party terms.
