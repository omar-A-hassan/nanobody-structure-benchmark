# 7R73 Docking Workflow

This directory contains the supplied input-preparation notebook for the illustrative 7R73 D7 VHH–HIV-1 gp120 docking case.

It prepares experimental and predicted VHH structures and standardizes the
inputs. It does not run HDOCK, DockQ, epitope scoring or the regional scorer.
It copies one existing IgFold PDB rather than executing the IgFold benchmark.

The runnable notebook has been cleared of execution counts and outputs,
without changes to its cell sources. The original executed copy remains in
Git history at the commit recorded in `provenance/docking_preparation_execution.json`.
Its install cells use a moving SimpleFold branch and unpinned dependencies;
this update does not claim an independently tested regeneration setup.
Prepared outputs are in `data/docking/`. Run regional rescoring locally with
`uv run --locked python src/score_prepared_vhh.py --check`.
