# Antigen and interface source overlap

This analysis compares antigen sequences, epitopes and heavy-domain contact
patterns between benchmark references and available model-linked experimental
source complexes. It measures source-complex similarity, not exact checkpoint
training membership.

## Run

From the repository root:

```sh
uv run --locked python src/audit_complex_overlap.py screen
uv run --locked python src/audit_complex_overlap.py contacts
```

The supplied coordinates permit contact replay without network access. To
rebuild the coordinate cache from the complete source archive, use:

```sh
uv run --locked python src/audit_complex_overlap.py contacts --archive ../sabdab2-splits.tar.gz
```

The archive is SAbDab2 ML v0.1.0, DOI 10.5281/zenodo.20083995. Its download
instruction is in the root README. The script verifies the recorded archive
size and SHA256 before reading exact member names. Hashed cache filenames
preserve case-sensitive source identities on case-insensitive filesystems.
Frozen benchmark coordinates and scores are not modified.

## Measurements

- Available source PDBs link IgFold and NanoBodyBuilder2 to source complexes.
  SimpleFold uses an antibody-complex proxy deposited by 1 May 2020, not a
  complete verified training list.
- Antigens use global affine-gap alignment, identical positions divided by
  all alignment columns, with a 90% identity threshold. Protein and peptide
  labels share an alphabet; DNA and RNA remain separate.
- Native contacts use heavy atoms within 5 angstrom, the first deposited model
  and positive occupancy. The highest-occupancy alternate atom is selected;
  source order resolves ties.
- Epitopes are antigen residues contacting the heavy variable domain. Interface
  patterns are the heavy-domain/antigen residue contact pairs, not the complete
  paired-antibody interface.
- Sequence alignments map residue indices between partners. Jaccard similarity
  includes unmapped positions in its union. One empty contact set gives zero
  overlap; two empty sets do not define a similarity.

## Results and coverage

`outputs/revision_2026-10-26/complex_overlap/summary.json` is the result summary.
The same directory supplies sequence matches, exact source-member identifiers,
source coordinates, native contacts, all comparison rows and per-target/model
summaries. `screen_inputs.json` records source hashes and thresholds.

The 588-target common set contains 281 retained polymer-antigen references;
280 support positive-occupancy native-contact measurement. One reference,
7N9TE, lacks positive-occupancy heavy-domain atoms. Other targets have only
nonpolymeric partners (91), no retained partner despite a bound annotation
(194), or a free annotation (22). They are not negative overlap results.

| Source pool | Antigen identity >=90% | Epitope Jaccard >=0.5 | Contact Jaccard >=0.5 |
|---|---:|---:|---:|
| SimpleFold proxy | 120 | 34 | 8 |
| IgFold | 132 | 88 | 40 |
| NanoBodyBuilder2 | 137 | 101 | 74 |

Counts are target-level maxima; one target can match several source complexes.
The 393-target VHH-screened subset has 192 polymer-antigen references and still
contains antigen matches. Thus VHH/H3 screening does not certify complex-level
independence. Per-target fields also separate same-PDB and different-PDB matches.
Nonpolymeric partners, source-removed partners and lower-identity structural
analogues are outside this sequence-based screen.

## Checks

```sh
uv run --locked python -m unittest discover -s tests -p test_complex_overlap.py -v
```

Regression checks cover alignment gaps, unmapped positions, empty contact sets,
hydrogens, first-model selection, zero occupancy and case-sensitive source
identity. The full repository suite has 12 tests. Source-contact availability
and source-chain failures are recorded explicitly rather than filled with zeros.
