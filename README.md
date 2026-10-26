# Nanobody structure benchmark

Code and data for the structural comparison of SimpleFold, IgFold and
NanoBodyBuilder2 on 597 experimental VHH references. The common comparison
contains 588 targets. The available-source training-overlap and redundancy
screen contains 393 representatives.

This package reproduces the supplied structural analyses. It does not yet
reproduce every result in the manuscript. The original IgFold inference
workflow and downstream docking evaluation still need integration. Prepared
docking inputs, a preparation notebook and regional rescoring are included. See
[Revision status](docs/REMAINING.md) lists completed responses and only three
remaining handover inputs: IgFold generation code/settings, docking evaluation
code/data, and the licence/DOI deposit.

## Reproduce the saved analyses

Install [uv](https://docs.astral.sh/uv/), open a terminal in this directory,
and check the supplied files:

```sh
uv run --locked python reproduce.py
```

Repeat the analyses and figures:

```sh
uv run --locked python reproduce.py --replay
```

The first command checks file hashes, target counts, reference reconstruction
and the retained prediction archives. The second also runs the statistical
analyses, overlap screens, reference checks, sampling figures and all 110 ESM2
probe fits. It compares numerical outputs with the supplied results.

Replay uses a temporary copy. It writes logs, comparison reports and generated
figures to `verification/`; it does not replace the primary scores or launch
model inference. `verification/reexecution.json` gives the result of each job.
PDF bytes can differ because of timestamps; numerical comparisons are checked
separately. Python 3.12 and 3.13 environments use the supplied lockfiles.
The first run may download Python and dependencies. Add `--offline` only after
these dependencies are cached. No system Python packages are installed.

Check the snapshot manifest before replay. A replay adds verification files;
keep the supplied manifest as the record of the distributed snapshot.

## Inputs and outputs

| Analysis | Inputs and code | Supplied results |
|---|---|---|
| Regional accuracy | `data/inputs/targets_vhh.json`, `outputs/structure/`, `src/reviewer_revision_analysis.py` | All-loop scores, distributions and paired intervals in `outputs/revision_2026-10-16/` |
| Training overlap | `src/revision_training_audit.py`, `src/revision_overlap_extension.py` and recorded source pools | Per-target matches, screened membership and summaries in `outputs/revision_2026-10-17/overlap_extension/` |
| Experimental references | `src/reconstruct_revision_references.py` and 604 bundled source coordinate members | Exact 597-record reconstruction and reference-policy checks |
| Sampling and probe | `outputs/{trajectory,hard,hitrate,ensemble}/`, `data/esm_probe/`, `src/reproduce_esm_probe.py` | Disagreement, finite-coverage and grouped probe results |
| Inference repeats | `kaggle/kernels/revision-archive-20261023/`, `inference/repeat_outputs/` | Complete SimpleFold v1 and NanoBodyBuilder2 v2 coordinates, scores, run records and comparisons |
| Prepared docking inputs | `data/docking/`, `src/score_prepared_vhh.py` | Four prepared-input H1/H2/H3/framework score rows and input hashes in `outputs/docking/` |
| Figures | The analysis drivers, `src/recreate_reused_figures.py` and `src/make_revision_workflow.py` | Regional, sampling, probe and workflow figures |

Existing analysis filenames and dated result paths are retained so scripts,
source records and result comparisons keep the same references.
[Dependencies](docs/DEPENDENCIES.md) lists each entry point and environment.
[Coverage](docs/COVERAGE.md) describes the supplied data and their limits.

## Reconstruct the complete source archive

The bundled coordinate members permit inspection without a full download.
To validate extraction from the complete SAbDab2 release:

```sh
curl --fail --location --output ../sabdab2-splits.tar.gz 'https://zenodo.org/records/20083995/files/splits.tar.gz?download=1'
uv run --locked python reproduce.py --replay --reference-archive ../sabdab2-splits.tar.gz
```

Replay verifies the recorded source checksum before extraction. Keep the
download outside the repository so it does not change the snapshot inventory.
The complete source archive is about 876 MB and is not included in this package.

## Run new predictions

The SimpleFold and NanoBodyBuilder2 repeat notebooks use the frozen
`data/inputs/targets_vhh.json`. For a new Kaggle run, upload this file as a
dataset, attach one copy, set the notebook metadata to the account and dataset,
enable internet, and select two T4 GPUs. The included metadata identifies the
completed runs; access to those private runs is not required to run a new copy.

Keep new coordinates, scores and run records in a separate version. The repeat
archives do not replace the primary benchmark. See
[Completed repeats](docs/COMPLETED_REPEATS.md) for score comparisons and local
coordinate-rescoring commands. These commands score retained coordinates;
they do not generate new structures.

The ESM2 probe replay fits classifiers from the supplied cached embeddings.
It does not run ESM2 embedding inference. The original sampling, trajectory
and hit-rate generators are included under `kaggle/kernels/`. Their saved
results support analysis replay; the recent replay did not rerun GPU sampling.
[Sampling](docs/SAMPLING.md) gives the notebook settings, inputs, execution
order and outputs. Full trajectory archives and new GPU sampling runs are not
required for this handover. Saved-data replay is not a new generation run.

## Check the prepared docking inputs

```sh
uv run --locked python src/score_prepared_vhh.py --check
```

This command verifies all four prepared VHH regional-score rows against the
frozen 7R73A reference. It reuses the original Kabsch function and checks the
fit independently with BioPython. Omit `--check` to regenerate the regional
CSV and JSON. It does not reproduce DockQ or epitope endpoints. The supplied
preparation notebook and docking summary CSVs are described in
[Docking case](docs/DOCKING.md). The missing downstream evaluation procedure
and its data are item D1 in [Revision status](docs/REMAINING.md).

## Tests and source records

```sh
uv run --locked --with-requirements environments/inference-check.lock python -m unittest discover -s tests -v
```

`MANIFEST.sha256` records file hashes. `provenance/source_selection.json`
records selected source artifacts and storage changes. Original numerical
outputs and original inference notebooks remain distinct from later repeats.
Pretrained model weights are not redistributed.

## Licence and citation

The project licence and public version DOI have not been assigned. Before
describing a deposit as the complete manuscript archive, add the remaining
[handover inputs](docs/REMAINING.md) and match its declared coverage to the
final paper.
Existing source-data DOIs identify those datasets, not this repository.
Upstream attribution and terms are retained in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
