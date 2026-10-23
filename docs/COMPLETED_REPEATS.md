# Completed monomer repeats

These are new predictions from the frozen 597-target VHH input file. The
primary reference and score collections are unchanged. All coordinate ZIP
members pass CRC, byte-count and SHA256 checks. The unchanged original scorer
reproduces every saved repeat row exactly in the local locked environment.

| Run | Coordinates | Whole-chain H3 median, common 588 | Median absolute H3 change | Maximum absolute H3 change |
|---|---:|---:|---:|---:|
| SimpleFold v1 | 597 | 2.944 Å | 0.000 Å | 0.001 Å |
| NanoBodyBuilder2 v1 | 588 | 2.474 Å | 0.018 Å | 1.967 Å |
| NanoBodyBuilder2 v2 | 588 | 2.474 Å | 0.014 Å | 1.952 Å |

Changes are paired against the primary scores. The primary common-set medians
are 2.944 Å for SimpleFold and 2.4835 Å for NanoBodyBuilder2. Differences between
medians and medians of within-target changes are not the same statistic.
Full eight-region/alignment summaries and per-target comparisons are supplied.
Small median changes do not imply identical structures: some NBB2 targets
change by almost 2 Å. Both NBB2 repeats return the same 588 IDs and refuse the
same nine incomplete inputs.

## Execution records

`inference/repeat_outputs/` contains the coordinate archives, source notebooks,
metadata, input hashes, environments, logs and separate comparisons. SimpleFold
v1 records source commit `466273494290b3dff842b8a78bdc861c2a642d6f`, its 100M
flow checkpoint hash, 500 steps, tau 0.01, seed 12345, one sample, no pLDDT,
two T4 workers, the shard enumeration and the 90-second startup stagger.
Its source TAR records the pinned commit before the explicitly recorded
one-line half-precision ESM patch. ESM2 is identified by model name; its
weight-file hash was not logged.

NBB2 v2 uses ImmuneBuilder 1.2 with default refinement, two T4 workers, installed
source hashes and all four numbered model-weight hashes. V1's empty weight-hash
map is preserved and explained separately. No weight files are redistributed.
Logged hashes are execution records, not independent rehashes of omitted weights.

Some inherited console summaries contain literature values or all-597
SimpleFold values beside 588 NBB2 outputs. Those console comparisons are not
the common-set manuscript analysis. Use the saved per-target records and
`src/summarize_completed_repeats.py`, which explicitly intersects all three
primary score sets to obtain 588 common targets.

## Local replay

```sh
uv run --no-project --python 3.13 --with-requirements environments/inference-check.lock python src/verify_archive_reruns.py --model simplefold --run-directory inference/repeat_outputs/simplefold_version1 --report verification/simplefold_version1_reexecution.json --require-exact --require-checkpoint-hashes
uv run --no-project --python 3.13 --with-requirements environments/inference-check.lock python src/verify_archive_reruns.py --model nbb2 --run-directory inference/repeat_outputs/nbb2_version2 --report verification/nbb2_version2_reexecution.json --require-exact --require-checkpoint-hashes
uv run --no-project python src/summarize_completed_repeats.py
```

These commands do not launch inference or contact Kaggle. Scoring uses a
disposable coordinate directory. The imported verification receipts are not
overwritten. Exact historical coordinate identity cannot be tested against
deleted coordinate files; numerical agreement is reported at saved precision.
The original IgFold inference package and docking inputs/poses are not supplied
by these monomer repeats.
