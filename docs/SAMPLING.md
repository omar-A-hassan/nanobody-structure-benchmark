# Sampling code and saved analyses

The repository includes generation notebooks, target lists, sampler settings
and saved sampling records. Local replay uses those records; generation is a
separate GPU workflow.

## Reproduce the reported analyses without inference

From the repository root:

```sh
uv run --locked python reproduce.py --replay
```

This uses the saved sampling records and cached ESM2 probe embeddings. It
recreates analyses and figures; it does not run folding or new ESM2 inference.

## Included generation notebooks

| Notebook under `kaggle/kernels/` | Generation settings | Main output in `/kaggle/working/` |
|---|---|---|
| `sf-traj-tfm-597/sf-traj-tfm-597.ipynb` | 597 targets, six samples each, 500 steps, tau 0.01, seed 12345; log every five steps and the final step | `traj_results.json`; per-target arrays in `traj_npy/` |
| `sf-hard/sf-hard.ipynb` | 24 hard and 12 easy targets, six samples each, 500 steps, tau 0.01, seed 12345 | `hard_results.json` |
| `sf-hitrate/sf-hitrate.ipynb` | 30 hard and ten easy targets, 20 samples each at tau 0.01, 0.8, 1.5 and 0.3, 500 steps, seed 12345 | `hitrate_results.json` |

For a new generation run, upload `data/inputs/targets_vhh.json` as a Kaggle
dataset, attach exactly one copy, enable internet, select two T4 GPUs and run
the selected notebook's cells in order through the main result-JSON writer.
The optional post-hoc TFM-model cell in the full trajectory notebook is not
needed for the reported spread/trajectory analysis. The notebooks install
SimpleFold and dependencies, retrieve model checkpoints and apply their
half-precision ESM patch. They use `simplefold_100M` and ESM2-3B. They sort
targets and alternate them between two workers, with a 60-second stagger.

The full trajectory notebook retains an earlier 48-target setup cell. Run
the later full-benchmark replacement cell as well; it selects all 597 inputs
before generation. Its logger saves per-residue means of unmasked atom
coordinates, not C-alpha-only positions.

The notebooks retain their recorded installation commands, including the
`fix-multichain-fasta-parsing` SimpleFold branch. A new run can install different
source revisions or dependency versions from those recorded in the original
execution. Keep new results separate from the supplied records. Notebook
Markdown descriptions have been clarified; code cells are unchanged.

The `sf-ensemble/` and `esm-probe/` folders contain the separate ensemble and
embedding-generation sources. The manuscript's analysis replay uses the
corresponding saved records, not those generation steps.
