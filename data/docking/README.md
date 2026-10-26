# 7R73 Docking Inputs

This directory contains the standardized structural inputs used in the illustrative 7R73 D7 VHH–HIV-1 gp120 re-docking case.

## Prepared inputs

- `7R73_antigen_DOCK.pdb` — standardized experimental gp120 antigen receptor
- `7R73_experimental_DOCK.pdb` — standardized experimental D7 VHH positive control
- `7R73A_IgFold_DOCK.pdb` — standardized IgFold-predicted VHH
- `7R73A_NBB2_DOCK.pdb` — standardized NanoBodyBuilder2-predicted VHH
- `7R73A_SimpleFold_DOCK.pdb` — standardized SimpleFold-predicted VHH

All four VHH docking inputs contain the same 126 residues and 973 heavy atoms after standardization. The antigen contains 287 resolved residues and 2283 heavy atoms.

Hydrogens were removed, protein ATOM records were retained, missing heavy atoms for residues already present were restored, unresolved missing residues were not rebuilt, the VHH chain identifier was normalized to H, and the antigen chain remained G.

The same antigen coordinates were used for all docking runs. HDOCK was run in template-free mode without receptor or ligand binding-site restraints.

The preparation notebook is in `kaggle/kernels/7r73-docking/`. Regional
rescoring against the frozen 7R73A reference is implemented in
`src/score_prepared_vhh.py`; the CSV/JSON in `outputs/docking/` records the
four VHH file hashes, atom counts and scores. The fixed antigen is covered
by the repository checksum manifest. The link from these supplied inputs to
the historical HDOCK submissions remains to be documented; see
`docs/REMAINING.md`.
