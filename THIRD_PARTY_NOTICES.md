# Third-party sources and transformations

The repository's original code and documentation are licensed under the MIT
License in `LICENSE`. Upstream datasets and third-party materials retain their
own source terms; this project licence does not relicense those materials or
clear separate model-weight terms. The notices below retain source attribution.

1. **SAbDab2 Machine Learning Dataset, v0.1.0.** Capel, Vavourakis, Williams,
   Taylor and Deane; University of Oxford (2026).
   Source record: <https://doi.org/10.5281/zenodo.20083995>.
   Retained official record metadata identifies CC BY 4.0:
   <https://creativecommons.org/licenses/by/4.0/>. Changes consist of benchmark
   filtering, sequence/span/coordinate extraction, alternate-atom comparison,
   statistical analyses and lossless compression of a source table. The 604
   coordinate members are a subset of the source archive. These transformations
   do not imply endorsement by the source authors. The upstream README and
   official record metadata are retained.

2. **Training data for IgFold, v1.** Ruffolo, Chu, Mahajan and Gray (2023).
   Source: <https://doi.org/10.5281/zenodo.7820263>. Official record metadata
   identifies CC BY 4.0. Changes comprise experimental heavy/H3 sequence-pool
   extraction and experimental/OAS sequence overlap screens. Original source
   archive checksums, OAS coverage/numbering evidence and retrieval receipts
   are retained; the source ZIPs are not redistributed here. This data licence
   does not describe the separate IgFold model-code or model-weight terms.
   One unmodified 7NXB coordinate member is included as an attributed H3-span
   regression fixture, with its archive/member receipt in `data/regression/`.

3. **ImmuneBuilder.** <https://github.com/oxpig/ImmuneBuilder>.
   BSD 3-Clause; Copyright (c) 2022 Brennan Abanades Kenyon. The complete
   upstream notice is retained in `provenance/ImmuneBuilder-LICENSE.txt`.
   The public `data/nanobody_data.csv` snapshot was retrieved on 16 October
   2026. It is evidence of the public training list, not verification of exact
   historical checkpoint membership. No model weights are redistributed.

4. **SimpleFold.** <https://github.com/apple-aiml-research/ml-simplefold>.
   Copyright (c) 2025 Apple Inc.; the full upstream permission notice is
   retained in `provenance/SimpleFold-LICENSE.txt`. Public recipe/manifests,
   repository metadata and evidence excerpts are retained with source URLs and
   receipts. The extended identifier dictionary is losslessly gzipped. Model
   and language-model weights are excluded and retain their separate terms.

5. **ANARCII 2.0.8.** <https://github.com/oxpig/ANARCII>.
   The complete retrieved BSD 3-Clause LICENCE, retrieval receipt and PyPI
   metadata are retained under
   `outputs/revision_2026-10-17/overlap_extension/sources/`. ANARCII numbering
   was used for candidate-domain/H3 recovery. No model-weight package is
   redistributed. The collection dependency lock is separate from analysis.

6. **RCSB/wwPDB.** Public entity records and obsolete mmCIF coordinates support
   source recovery. <https://www.wwpdb.org/about/usage-policies> identifies
   PDB archive data as CC0 1.0. Source policies, record URLs, retrieval receipts
   and sequence/coordinate derivations are retained in the overlap cache.

7. **UniProt.** <https://www.uniprot.org/help/license>.
   Copyrightable database content is attributed under CC BY 4.0 in the retained
   official policy evidence. Three current canonical FASTAs and source receipts
   support substring/global-alignment checks of entity-level accession flags.
   These checks do not validate historical AFDB model sequences.

The listed DOIs identify upstream datasets. None is a DOI for this repository.
Public source receipts describe the verified retrieval context and are not a
claim of complete exposure coverage or author endorsement.
