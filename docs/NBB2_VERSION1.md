# Completed NanoBodyBuilder2 repeat, version 1

The frozen run is in `inference/repeat_outputs/nbb2_version1/`. It was copied
from the study's `outputs/revision_2026-10-23/repeat_runs/nbb2_version1/`, not
from version 1 rather than the later version on the same Kaggle reference,
`potatoguy44/nbb2-benchmark-archive-20261023`.

## Verified facts

| Check | Result |
|---|---|
| Run status receipt | COMPLETE |
| Retained coordinates | 588 PDB members in `output/nbb2_coordinates.zip` |
| Coordinate ZIP SHA256 | `e7d0c8e1c586e4372967d1dcca4041ca2687110e4e2190533705cdf1bf65e93e` |
| ZIP CRC, manifest names/sizes/SHA256 | Passed for every member |
| Primary/repeat target IDs | All 588 identical; no original-only or repeat-only IDs |
| Saved coordinate-to-score replay | Complete rows exactly equal; maximum regional-score difference 0.0 |
| Whole-chain H3 median RMSD | Primary 2.4835 Å; repeat 2.474 Å |
| Whole-chain H3 median absolute per-target change | 0.018 Å |
| Whole-chain H3 maximum absolute per-target change | 1.967 Å |
| Primary benchmark replacement | No |
| V1 checkpoint hashes | Not recorded; raw map is `{}` |

The separate `comparison.json` contains all eight regional/alignment summaries;
`comparison.csv` contains all 4,704 per-target regional comparisons. The H3
figures above refer to whole-chain alignment. Framework-aligned H3 maximum
change is 2.052 Å. Differences between medians and medians of paired changes
are distinct statistics.

There are 597 frozen input references; the nine absent NBB2 targets are
10ZOA, 7A5VO, 9FAKG, 7QHAC, 7PQQB2, 8AFEJ, 8XUMC, 8JJVA and 9PYMB,
matching the nine absent primary NBB2 scores. Scoring the retained 588 succeeds.
The new coordinates are repeat predictions, not recovered historical
predictions. Exact rescoring of repeat coordinates does not establish identity
with the original deleted prediction coordinates.

## Incomplete checkpoint provenance

V1's notebook only selected checkpoint paths ending in `.pt` or `.pth`.
ImmuneBuilder's `nanobody_model_*` files are extensionless, so the retained
`output/run_provenance.json` has `checkpoint_sha256: {}`. This is a logging gap,
not a coordinate or scoring failure. No hashes are backfilled from the earlier
three-target CPU check, a later checkpoint download, or the running v2.

The archived v1 `repeat.ipynb` preserves that original filter. The separately
included current `src/prepare_archive_reruns.py` fixes the filter and asserts
that all four numbered NBB2 weights are logged. That source was not run during
the v1 package extension. The subsequently completed NBB2 v2 and SimpleFold v1
are retained separately; see `COMPLETED_REPEATS.md`. They do not replace or
backfill v1 artifacts.

V1's environment receipt is preserved rather than replaced by the local check
environment. It records Python 3.13.15, ImmuneBuilder 1.2, OpenMM 8.6.1 and
Biotite 2.0.0. Local rescoring used the separately locked NumPy 2.5.4,
SciPy 1.18.1 and Biotite 1.8.0 environment and reproduced the saved rows exactly.

## Local recheck without inference or cloud requests

```sh
uv run --no-project --python 3.13 --with-requirements environments/inference-check.lock python src/verify_archive_reruns.py --model nbb2 --run-directory inference/repeat_outputs/nbb2_version1 --report verification/nbb2_version1_reexecution.json --require-exact
```

The original scorer is extracted from the recovered NBB2 benchmark notebook.
The verifier checks the frozen input hash and every ZIP member before scoring,
uses a disposable extraction directory, and writes the new report separately.
The imported `inference/repeat_outputs/nbb2_version1/verification.json` remains
the original saved verification. No 588-file uncompressed duplicate is shipped.

For later runs expected to record weights, add `--require-checkpoint-hashes`.
It asserts complete SHA256 entries for all four numbered NBB2 weights; the flag
must be omitted for v1. Logged weight hashes are not independently rehashed
without the weight files, which are excluded from this repository.

The verifier's optional `--download` mode pulls source plus metadata and outputs
after the remote state is COMPLETE. It addresses the current remote version,
so it cannot be used to refresh this frozen v1 snapshot. The local versioned
mode excludes `--download` and preserves the original verification receipt.

The preceding 15-check `verification/reexecution.json` and its associated
artifacts are unchanged. The extension receipt is
`provenance/nbb2_version1_extension.json`; the preceding public checksum context
is retained under `provenance/core_snapshot/`.
