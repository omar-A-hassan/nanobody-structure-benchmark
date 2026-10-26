# Illustrative 7R73 Docking Case

## Overview

The downstream analysis uses the experimentally determined D7 VHH–HIV-1 gp120 complex from PDB 7R73 as an illustrative controlled re-docking case.

The analysis compares four VHH inputs:

1. Experimental D7 VHH
2. NanoBodyBuilder2 prediction
3. IgFold prediction
4. SimpleFold prediction

The experimentally observed gp120 conformation was held fixed across all docking runs.

This is a single-complex illustrative case and is not intended to establish a general relationship between VHH regional prediction error and docking performance.

The supplied files support input inspection and regional monomer rescoring.
The HDOCK/DockQ protocol and endpoint summaries below are reported records;
they have not been reproduced locally from retained docking poses. Missing
pose-scoring material is listed in [Remaining work](REMAINING.md).

## Case Selection

7R73 was selected as an informative pilot because it provides:

- an experimentally resolved VHH–antigen complex,
- a well-defined native interface,
- a high-resolution X-ray structure at 1.76 Å,
- sufficient resolved VHH and antigen coordinates for controlled re-docking and residue-level epitope evaluation,
- and a substantial contrast in predicted CDR-H3 accuracy among the evaluated models.

The case contains the D7 nanobody bound to HIV-1 gp120.

## Prepared Structural Inputs

The exact standardized docking inputs are stored in:

`data/docking/`

Files:

- `7R73_antigen_DOCK.pdb`
- `7R73_experimental_DOCK.pdb`
- `7R73A_IgFold_DOCK.pdb`
- `7R73A_NBB2_DOCK.pdb`
- `7R73A_SimpleFold_DOCK.pdb`

All four VHH docking inputs were standardized to contain:

- 126 residues
- 973 heavy atoms
- chain identifier H

The antigen contains:

- 287 resolved residues
- 2283 heavy atoms
- chain identifier G

During preparation:

- hydrogens were removed,
- protein ATOM records were retained,
- missing heavy atoms belonging to existing residues were restored,
- unresolved missing residues were not reconstructed,
- VHH chain identifiers were normalized to H,
- the antigen chain remained G.

## Docking Protocol

Docking was performed using the HDOCK web server in template-free mode.

The same standardized gp120 receptor structure was used for every VHH input.

No receptor or ligand binding-site restraints were supplied.

For each VHH input, the 10 highest-ranked poses according to HDOCK were retained.

The experimental VHH served as a positive-control re-docking input.

Because the experimentally observed antigen conformation was retained, this experiment is a holo re-docking case rather than fully unbound docking.

## Docking Evaluation

Docked poses were compared with the experimental 7R73 complex using DockQ v2.1.3.

The reported docking metrics include:

- DockQ
- ligand RMSD (LRMSD)
- interface RMSD (iRMSD)
- fraction of native contacts (fnat)

Two pose-selection summaries were used:

### Top-1

The pose ranked first by HDOCK.

### Best-of-Top-10

The pose with the highest DockQ among the 10 retained HDOCK candidates.

Because the experimental structure is required to identify the highest-DockQ candidate, Best-of-Top-10 is treated as an oracle measure of sampling capacity rather than a deployable ranking method.

## Epitope Recovery

The experimental gp120 epitope was defined as antigen residues containing at least one heavy atom within 5 Å of any VHH heavy atom in the native complex.

Predicted epitopes were generated from docked poses using the same 5 Å heavy-atom cutoff.

Predicted and experimental epitope sets were compared using:

- recall
- precision
- Jaccard similarity

The experimental gp120 epitope contains 26 residues.

## Main Top-1 Results

| VHH input | DockQ | LRMSD (Å) | fnat | Epitope recovered | Recall | Jaccard |
|---|---:|---:|---:|---:|---:|---:|
| Experimental | 0.986 | 0.640 | 1.000 | 26/26 | 1.000 | 0.963 |
| NanoBodyBuilder2 | 0.740 | 2.926 | 0.700 | 20/26 | 0.769 | 0.645 |
| IgFold | 0.022 | 35.767 | 0.000 | 2/26 | 0.077 | 0.042 |
| SimpleFold | 0.009 | 60.919 | 0.000 | 0/26 | 0.000 | 0.000 |

The supplied top-1 and best-of-top-10 summary tables are stored in:

`outputs/docking/`

## Regional Structural Accuracy

Framework-aligned CDR-H1, CDR-H2, CDR-H3, and framework RMSDs are reported for the exact docking VHH inputs using the structural scoring procedure described in the manuscript.

The supplied prepared files give:

| VHH input | H1 (Å) | H2 (Å) | H3 (Å) | Framework (Å) |
|---|---:|---:|---:|---:|
| NanoBodyBuilder2 | 2.660 | 0.724 | 1.059 | 0.928 |
| IgFold | 3.079 | 1.424 | 2.613 | 0.889 |
| SimpleFold | 2.196 | 0.813 | 5.763 | 0.960 |
| Experimental control | 0.000 | 0.000 | 0.000 | 0.000 |

Run `uv run --locked python src/score_prepared_vhh.py --check` from the
repository root. The scorer verifies identical 126-residue sequences, uses
492 matched N/CA/C/CB atoms and fits 350 framework atoms to the frozen 7R73A
reference. No CDR is independently fitted. The JSON and CSV retain input
hashes, full-precision regional scores and an independent BioPython alignment
check. The control's zero values are computed monomer RMSDs, not a claim of
perfect re-docking. These scores do not establish which files were submitted
to HDOCK; the submission linkage is a separate provenance task.

`verification/prepared_vhh_validation.json` records a fresh distributed-file
copy check: package inventory, regional-score replay and all eight unit tests.
It is separate from the earlier full saved-analysis replay and from inference.

## Notebook

The docking preparation workflow is stored under:

`kaggle/kernels/7r73-docking/`

The notebook generates and standardizes the structural inputs. It copies an
existing IgFold prediction and generates SimpleFold/NanoBodyBuilder2 inputs.
It does not run the 597-target IgFold benchmark, HDOCK, DockQ, regional scoring
or epitope scoring, and it does not generate the supplied docking summary CSVs.
The runnable copy has no saved outputs. Its executed source is retained at
the commit identified in `provenance/docking_preparation_execution.json`.

## Scope and Limitations

This analysis contains one nanobody–antigen complex.

Accordingly:

- no cross-complex statistical relationship is inferred,
- the analysis does not isolate the effect of any single CDR,
- docking performance may reflect regional geometry, framework geometry, structure preparation, docking sampling, and ranking,
- the antigen was retained in its experimentally observed bound conformation,
- HDOCK was the only docking method evaluated,
- no between-submission variance or independent docking/scoring control was estimated,
- docking and epitope agreement do not establish binding affinity, neutralization, or biological activity.

The results should therefore be interpreted as an illustrative case study.
