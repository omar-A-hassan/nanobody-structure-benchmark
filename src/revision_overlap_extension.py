"""Source-bound B2 extension; never imports the side-effectful 16 October audit.

Run with isolated uv dependencies. All generated data stays in overlap_extension.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.request import Request, urlopen
from urllib.parse import urlencode
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import re
import time
import zipfile
import gzip
import io
from collections import Counter
import html

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "outputs/revision_2026-10-16"
OUT = ROOT / "outputs/revision_2026-10-17/overlap_extension"
SRC = OUT / "sources"
csv.field_size_limit(10_000_000)


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def save_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n")


def save_csv(path, rows):
    assert rows
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def frozen_inputs():
    """Verify the explicitly curated PUBLIC subset, including compressed sources."""
    receipt = OUT / "frozen_input_sha256.json"
    expected = json.loads(receipt.read_text())
    now = {name: digest(ROOT / name) for name in expected}
    assert now == expected, "Public frozen inputs changed; use a fresh checkout for replay"
    return now


def fetch(url, name, md5=None, size=None):
    """Cache complete validated files; retain interrupted downloads for resumption."""
    path = SRC / name
    receipt = path.with_name(path.name + ".receipt.json")
    if path.exists():
        assert size is None or path.stat().st_size == size
        assert md5 is None or digest(path, "md5") == md5
        assert receipt.exists() and json.loads(receipt.read_text())["sha256"] == digest(path)
        print("cached", name, path.stat().st_size, flush=True)
        return path
    partial = path.with_name(path.name + ".part")
    offset = partial.stat().st_size if partial.exists() else 0
    headers = {"User-Agent": "B2-public-source-overlap-audit/1.0"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    started = time.monotonic()
    with urlopen(Request(url, headers=headers), timeout=120) as response:
        append = offset > 0 and response.status == 206
        if append:
            assert response.headers["Content-Range"].startswith(f"bytes {offset}-")
        with partial.open("ab" if append else "wb") as f:
            last = time.monotonic()
            for chunk in iter(lambda: response.read(1024 * 1024), b""):
                f.write(chunk)
                if time.monotonic() - last > 20:
                    print("download", name, f.tell(), "bytes", flush=True)
                    last = time.monotonic()
        info = {"url": url, "resolved_url": response.url,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "content_type": response.headers.get("Content-Type"),
                "etag": response.headers.get("ETag")}
    assert partial.stat().st_size > 0
    assert size is None or partial.stat().st_size == size
    assert md5 is None or digest(partial, "md5") == md5, f"Checksum mismatch: {name}"
    partial.replace(path)
    info.update(bytes=path.stat().st_size, sha256=digest(path), expected_md5=md5,
                expected_bytes=size, elapsed_seconds=round(time.monotonic() - started, 2))
    save_json(receipt, info)
    print("fetched", name, info["bytes"], flush=True)
    return path


def collect_oas():
    meta = json.loads(fetch("https://zenodo.org/api/records/7820263", "zenodo7820263.json").read_text())
    sabdab_meta = json.loads(fetch("https://zenodo.org/api/records/20083995", "zenodo20083995.json").read_text())
    save_json(OUT / "zenodo_licenses.json", {str(m["id"]): {
        "doi": m.get("doi"), "license": m["metadata"].get("license"),
        "access_right": m["metadata"].get("access_right"),
        "source": f"https://zenodo.org/api/records/{m['id']}"} for m in (meta, sabdab_meta)})
    files = [f for f in meta["files"] if f["key"].startswith("af_oas_")]
    assert len(files) == 2
    def get(f):
        return fetch(f["links"]["self"], f["key"], f["checksum"].split(":")[1], f["size"])
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(get, files))


def collect_rcsb(benchmark=False):
    baseline = json.loads((OLD / "training_overlap_audit.json").read_text())
    missing = (sorted({t["pdb"].upper() for t in json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())})
               if benchmark else baseline["nbb_source"]["unmapped_pdbs"])
    if not benchmark:
        assert len(missing) == 247
    prefix = "benchmark" if benchmark else "rcsb"
    entities = []
    absent = []
    for start in range(0, len(missing), 16):
        batch = missing[start:start + 16]
        query = "{entries(entry_ids:" + json.dumps(batch) + "){rcsb_id polymer_entities{" + (
            "rcsb_id rcsb_polymer_entity{pdbx_description} "
            "rcsb_polymer_entity_container_identifiers{auth_asym_ids asym_ids reference_sequence_identifiers{database_name database_accession}} "
            "entity_poly{type pdbx_seq_one_letter_code_can}}}}")
        path = fetch("https://data.rcsb.org/graphql?" + urlencode({"query": query}), f"{prefix}_batch_{start:03}.json")
        response = json.loads(path.read_text())
        assert not response.get("errors"), response.get("errors")
        found = {e["rcsb_id"] for e in response["data"]["entries"] if e}
        absent += sorted(set(batch) - found)
        for entry in response["data"]["entries"]:
            if entry:
                entities.extend(dict(pdb=entry["rcsb_id"], **e) for e in entry["polymer_entities"] or [])
        print(prefix, "entries", start + len(batch), "/", len(missing), flush=True)
    save_json(OUT / ("benchmark_rcsb_entities.json" if benchmark else "nbb2_rcsb_entities.json"),
              {"requested_pdbs": missing, "absent_pdbs": absent, "entities": entities})


def collect_obsolete():
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    data = json.loads((OUT / "nbb2_rcsb_entities.json").read_text())
    entities, failures = [], []
    for pdb in data["absent_pdbs"]:
        url = f"https://files.wwpdb.org/pub/pdb/data/structures/obsolete/mmCIF/{pdb[1:3].lower()}/{pdb.lower()}.cif.gz"
        try:
            path = fetch(url, f"obsolete_{pdb}.cif.gz")
            d = MMCIF2Dict(io.StringIO(gzip.decompress(path.read_bytes()).decode()))
            descriptions = dict(zip(d["_entity.id"], d["_entity.pdbx_description"]))
            for i, eid in enumerate(d["_entity_poly.entity_id"]):
                entities.append({"pdb": pdb, "rcsb_id": f"{pdb}_{eid}", "obsolete": True,
                    "rcsb_polymer_entity": {"pdbx_description": descriptions[eid]},
                    "rcsb_polymer_entity_container_identifiers": {
                        "auth_asym_ids": d["_entity_poly.pdbx_strand_id"][i].split(","), "asym_ids": [],
                        "reference_sequence_identifiers": []},
                    "entity_poly": {"type": d["_entity_poly.type"][i],
                        "pdbx_seq_one_letter_code_can": "".join(d["_entity_poly.pdbx_seq_one_letter_code_can"][i].split())}})
        except Exception as e:
            failures.append({"pdb": pdb, "url": url, "error": str(e)})
            print("obsolete unavailable", pdb, str(e), flush=True)
    save_json(OUT / "nbb2_obsolete_entities.json", {"entities": entities, "failures": failures})


def collect_simplefold():
    sources = {
        "simplefold_README.md": "https://raw.githubusercontent.com/apple/ml-simplefold/main/README.md",
        "simplefold_tree.json": "https://api.github.com/repos/apple/ml-simplefold/git/trees/main?recursive=1",
        "simplefold_paper.html": "https://arxiv.org/html/2509.18480v1",
        "simplefold_pdb_sp_afesme.yaml": "https://raw.githubusercontent.com/apple/ml-simplefold/main/configs/data/pdb_sp_afesme.yaml",
        "simplefold_train.yaml": "https://raw.githubusercontent.com/apple/ml-simplefold/main/configs/experiment/train.yaml",
        **{name: "https://ml-site.cdn-apple.com/models/simplefold/" + name for name in (
            "swissprot_list.csv", "afesm_list.csv", "afesme_dict.json")}}
    def get(item):
        name, url = item
        try:
            fetch(url, name)
            return {"name": name, "url": url, "status": "retrieved"}
        except Exception as e:
            print("source unavailable", name, type(e).__name__, str(e), flush=True)
            return {"name": name, "url": url, "status": "unavailable", "error": str(e)}
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(get, sources.items()))
    save_json(OUT / "simplefold_collection.json", results)


def imgt_domain(result):
    numbering = result.get("numbering") or []
    vh = "".join(aa for _, aa in numbering if aa != "-")
    h3 = "".join(aa for (n, _), aa in numbering if 105 <= int(n) <= 117 and aa != "-")
    accepted = (result.get("chain_type") == "H" and result.get("score", 0) >= 19
                and 70 <= len(vh) <= 180 and bool(h3))
    return vh, h3, accepted


def number_nbb():
    from anarcii import Anarcii
    import importlib.metadata
    data = json.loads((OUT / "nbb2_rcsb_entities.json").read_text())
    obsolete = json.loads((OUT / "nbb2_obsolete_entities.json").read_text())
    entities = data["entities"] + obsolete["entities"]
    proteins = {e["rcsb_id"]: "".join(e["entity_poly"]["pdbx_seq_one_letter_code_can"].split())
                for e in entities if e["entity_poly"]["type"] == "polypeptide(L)"}
    # Number every protein entity, including antigen/fusion candidates, to expose ambiguity.
    numbered_path = OUT / "nbb2_anarcii_numbering.json"
    model = Anarcii(cpu=True, ncpu=4, batch_size=16)
    if numbered_path.exists():
        numbered = json.loads(numbered_path.read_text())
    else:
        numbered = {}
        items = list(proteins.items())
        for start in range(0, len(items), 48):
            numbered.update(model.number(dict(items[start:start + 48]), scfv=True))
            print("numbered NBB2 entities", min(start + 48, len(items)), "/", len(items), flush=True)
        save_json(numbered_path, numbered)
    candidates = []
    by_id = {e["rcsb_id"]: e for e in entities}
    for name, result in numbered.items():
        eid = name if name in by_id else name.rsplit("-", 1)[0]
        assert eid in by_id
        vh, h3, accepted = imgt_domain(result)
        e = by_id[eid]
        candidates.append({"pdb": e["pdb"], "entity_id": eid, "domain_id": name,
            "auth_chains": ",".join(e["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"] or []),
            "description": e["rcsb_polymer_entity"]["pdbx_description"], "obsolete": e.get("obsolete", False),
            "entity_sequence": proteins[eid], "vhh": vh, "cdr3": h3,
            "chain_type": result.get("chain_type"), "numbering_score": result.get("score"),
            "accepted_heavy_candidate": bool(accepted), "error": result.get("error"),
            "membership_uncertainty": "PDB-level manifest only; heavy-domain inference does not verify VHH identity or historical training chain"})
    save_csv(OUT / "nbb2_entity_domains.csv", candidates)
    accepted = [r for r in candidates if r["accepted_heavy_candidate"]]
    pdblist = []
    for pdb in data["requested_pdbs"]:
        es = [e for e in entities if e["pdb"] == pdb]
        cs = [r for r in accepted if r["pdb"] == pdb]
        pdblist.append({"pdb": pdb, "n_polymer_entities": len(es), "n_heavy_domains": len(cs),
            "candidate_entity_ids": ";".join(sorted({r["entity_id"] for r in cs})),
            "candidate_author_chains": ";".join(r["auth_chains"] for r in cs),
            "obsolete_archive_used": any(e.get("obsolete", False) for e in es),
            "status": "heavy_candidates_not_confirmed_training_chains" if cs else "no_accepted_heavy_domain"})
    save_csv(OUT / "nbb2_pdb_coverage.csv", pdblist)
    targets = json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())
    target_path = OUT / "benchmark_anarcii_numbering.json"
    if target_path.exists():
        target_numbers = json.loads(target_path.read_text())
    else:
        target_numbers = model.number({t["pdb"] + t["H"]: t["hseq"] for t in targets})
        save_json(target_path, target_numbers)
    validation = []
    for t in targets:
        tid = t["pdb"] + t["H"]
        vh, h3, ok = imgt_domain(target_numbers[tid])
        validation.append({"tid": tid, "heavy_candidate_accepted": bool(ok),
            "h3_matches_frozen_imgt": h3 == t["CDRH3"], "vh_matches_frozen": vh == t["hseq"],
            "frozen_cdr3": t["CDRH3"], "numbered_cdr3": h3})
    save_csv(OUT / "numbering_validation.csv", validation)
    summary = {"anarcii_version": importlib.metadata.version("anarcii"), "protein_entities_numbered": len(proteins),
        "pdbs_entity_sequences_recovered": len({e["pdb"] for e in entities}),
        "new_pdbs_with_heavy_candidates": len({r["pdb"] for r in accepted}),
        "accepted_heavy_domains": len(accepted), "unique_vh_h3_pairs": len({(r["vhh"], r["cdr3"]) for r in accepted}),
        "no_accepted_heavy_domain": [r["pdb"] for r in pdblist if not r["n_heavy_domains"]],
        "benchmark_h3_concordant": sum(r["h3_matches_frozen_imgt"] for r in validation),
        "benchmark_vh_concordant": sum(r["vh_matches_frozen"] for r in validation)}
    save_json(OUT / "nbb2_recovery_summary.json", summary)
    print(json.dumps(summary), flush=True)


def make_aligner():
    from Bio.Align import PairwiseAligner
    return PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                          open_gap_score=-5, extend_gap_score=-.5)


def identity(aligner, a, b):
    if not a or not b:
        return 0.
    counts = aligner.align(a, b)[0].counts()
    return counts.identities / (counts.identities + counts.mismatches + counts.gaps)


def init_nearest(pool):
    import numpy as np
    global ALIGNER, POOLS
    ALIGNER = make_aligner()
    POOLS = []
    for field in ("vhh", "cdr3"):
        unique = {}
        for r in pool:
            if r[field]:
                unique.setdefault(r[field], r["source_id"])
        seqs = sorted(unique)
        counts = np.array([[s.count(chr(i)) for i in range(65, 91)] for s in seqs])
        lengths = np.array([len(s) for s in seqs])
        POOLS.append((seqs, [unique[s] for s in seqs], counts, lengths))


def nearest_one(target):
    import numpy as np
    row = {"tid": target["pdb"] + target["H"]}
    for field, sequence in (("vhh", target["hseq"]), ("cdr3", target["CDRH3"])):
        seqs, names, counts, lengths = POOLS[0 if field == "vhh" else 1]
        histogram = np.array([sequence.count(chr(i)) for i in range(65, 91)])
        # Any alignment has <= histogram-intersection matches and >= max(lengths) columns.
        upper = np.minimum(counts, histogram).sum(axis=1) / np.maximum(lengths, len(sequence))
        order = np.argsort(-upper, kind="stable")
        best, best_id, n_align = 0., "", 0
        for ix in order:
            if upper[ix] <= best:
                break
            value = identity(ALIGNER, sequence, seqs[ix])
            n_align += 1
            if value > best:
                best, best_id = value, names[ix]
        row.update({f"max_{field}_identity": best, f"nearest_{field}_source": best_id,
                    f"n_{field}_alignments": n_align})
    return row


def compare_pool(label, pool):
    from concurrent.futures import ProcessPoolExecutor
    pool = sorted(pool, key=lambda r: r["source_id"])
    path = OUT / f"{label}_sequence_pool.csv"
    save_csv(path, pool)
    target_path = ROOT / "data/inputs/targets_vhh.json"
    signature = {"pool_sha256": digest(path), "targets_sha256": digest(target_path),
                 "method": "Bio1.87_global_2_-1_-5_-.5_histogram_bound_v1"}
    cache = OUT / f"{label}_nearest.csv"
    receipt = OUT / f"{label}_nearest.receipt.json"
    if cache.exists() and receipt.exists() and json.loads(receipt.read_text()) == signature:
        print("cached comparisons", label, flush=True)
        return read_csv(cache)
    targets = json.loads(target_path.read_text())
    rows = []
    with ProcessPoolExecutor(max_workers=4, initializer=init_nearest, initargs=(pool,)) as executor:
        for i, row in enumerate(executor.map(nearest_one, targets, chunksize=8), 1):
            rows.append(row)
            if i % 48 == 0:
                print("comparison", label, i, "/597", flush=True)
    assert len(rows) == 597 and len({r["tid"] for r in rows}) == 597
    save_csv(cache, rows)
    save_json(receipt, signature)
    return rows


def compare_nbb():
    domains = read_csv(OUT / "nbb2_entity_domains.csv")
    pool = [{"source_id": r["domain_id"], "vhh": r["vhh"], "cdr3": r["cdr3"],
             "source": "RCSB_or_obsolete_wwPDB_entity_ANARCII_heavy_candidate"}
            for r in domains if r["accepted_heavy_candidate"] == "True"]
    compare_pool("nbb2_extension", pool)


def parse_heavy_fasta(text):
    heavy, active = "", False
    for line in text.splitlines():
        if line.startswith(">"):
            active = line.rsplit(":", 1)[-1] == "H"
        elif active:
            heavy += line.strip()
    return heavy


def pdb_h3(text):
    from Bio.Data.PDBData import protein_letters_3to1_extended
    seen, h3 = set(), ""
    for line in text.splitlines():
        if line.startswith("ATOM") and line[21] == "H" and line[12:16].strip() == "CA":
            key = line[22:27]
            if key not in seen and 93 <= int(line[22:26]) <= 102:
                h3 += protein_letters_3to1_extended[line[17:20]]
                seen.add(key)
    return h3


def compare_oas():
    for label in ("af_oas_paired", "af_oas_unpaired"):
        path = SRC / (label + ".zip")
        assert path.exists(), f"Await validated download: {path}"
        pool, missing, resource_forks, fallback, light_only = [], [], [], [], []
        n_fasta = 0
        with zipfile.ZipFile(path) as archive:
            names = set(archive.namelist())
            for name in sorted(names):
                if not name.endswith(".fasta"):
                    continue
                if "__MACOSX" in Path(name).parts or Path(name).name.startswith("._"):
                    resource_forks.append(name)
                    continue
                n_fasta += 1
                fasta = archive.read(name).decode()
                heavy = parse_heavy_fasta(fasta)
                if not heavy:
                    suffixes = {line.rsplit(":", 1)[-1] for line in fasta.splitlines() if line.startswith(">")}
                    assert suffixes == {"L"}, f"Unclassified FASTA cannot be silently excluded: {name}"
                    light_only.append(name)
                    continue
                pdb_name = str(Path(name).with_suffix(".pdb"))
                if pdb_name not in names:
                    missing.append(name)
                    h3, coordinate_status = "", "missing_coordinate_file"
                else:
                    coordinate_text = archive.read(pdb_name).decode()
                    coordinate_status = "html_numbering_error" if "<html" in coordinate_text.lower() else "pdb"
                    h3 = "" if coordinate_status == "html_numbering_error" else pdb_h3(coordinate_text)
                if h3 and h3 not in heavy:
                    h3, coordinate_status = "", "coordinate_h3_not_in_fasta"
                pool.append({"source_id": name, "vhh": heavy, "cdr3": h3,
                    "source": "IgFold_Zenodo7820263_" + label,
                    "coordinate_status": coordinate_status,
                    "h3_method": "coordinate_Chothia93-102_IMGT-equivalent" if h3 else "unavailable"})
                if not h3:
                    fallback.append(len(pool) - 1)
                if len(pool) % 2000 == 0:
                    print("extracted", label, len(pool), flush=True)
        # Some published .pdb files contain AbNum error HTML; retain their FASTA
        # in the conservative sequence screen and recover H3 by explicit numbering.
        if fallback:
            from anarcii import Anarcii
            cache = OUT / f"{label}_fallback_numbering.json"
            if cache.exists():
                numbered = json.loads(cache.read_text())
            else:
                model = Anarcii(cpu=True, ncpu=4, batch_size=16)
                items = [(pool[i]["source_id"], pool[i]["vhh"]) for i in fallback]
                numbered = {}
                for start in range(0, len(items), 48):
                    numbered.update(model.number(dict(items[start:start + 48])))
                    print("numbered OAS fallback", label, min(start + 48, len(items)), "/", len(items), flush=True)
                save_json(cache, numbered)
            for i in fallback:
                row = pool[i]
                _, h3, accepted = imgt_domain(numbered[row["source_id"]])
                if accepted:
                    row["cdr3"], row["h3_method"] = h3, "FASTA_ANARCII_2.0.8_IMGT105-117"
        assert n_fasta == len(pool) + len(light_only)
        save_json(OUT / f"{label}_coverage.json", {"archive": str(path.relative_to(ROOT)),
            "fasta_files": n_fasta, "heavy_pairs": len(pool), "unique_heavy_h3_pairs": len({(r["vhh"], r["cdr3"]) for r in pool}),
            "light_only_fasta_members": len(light_only),
            "missing_coordinate_pairs": missing, "archive_sha256": digest(path),
            "ignored_resource_fork_fasta_members": len(resource_forks),
            "coordinate_status_counts": dict(Counter(r["coordinate_status"] for r in pool)),
            "h3_method_counts": dict(Counter(r["h3_method"] for r in pool)),
            "unavailable_h3_members": [r["source_id"] for r in pool if not r["cdr3"]],
            "h3_definition": "Chothia93-102 with insertion positions; same span as 16 October IgFold audit"})
        compare_pool(label, pool)


def self_test():
    import random
    aligner = make_aligner()
    assert identity(aligner, "AAAA", "AAA") == .75
    assert identity(aligner, "ACDEFG", "ACDEFG") == 1
    rng = random.Random(20261017)
    pool = [{"source_id": str(i), "vhh": "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=rng.randint(5, 35))),
             "cdr3": "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=rng.randint(2, 12)))} for i in range(100)]
    init_nearest(pool)
    for i in range(20):
        target = {"pdb": str(i), "H": "H", "hseq": pool[i]["vhh"][:-1] + "A", "CDRH3": pool[i]["cdr3"] + "G"}
        fast = nearest_one(target)
        for key, query in (("vhh", target["hseq"]), ("cdr3", target["CDRH3"])):
            assert fast[f"max_{key}_identity"] == max(identity(aligner, query, r[key]) for r in pool)
    fixture = ROOT / "data/regression/igfold_7nxb.pdb"
    receipt = json.loads(fixture.with_suffix(".receipt.json").read_text())
    assert digest(fixture) == receipt["sha256"]
    assert pdb_h3(fixture.read_text()) == "AKDGGKLWVYYFDY"
    assert parse_heavy_fasta(">a:L\nGG\n>a:H\nAC\nDE\n") == "ACDE"
    assert parse_heavy_fasta(">:L\nACDE\n") == ""
    assert imgt_domain({"numbering": [], "chain_type": None}) == ("", "", False)
    engineered = {"numbering": [((n, " "), "A") for n in range(1, 129)], "chain_type": "H", "score": 29}
    assert imgt_domain(engineered)[2], "High-score engineered domains need not conserve canonical anchors"
    save_json(OUT / "self_test.json", {"status": "passed", "checks": ["global_alignment_definition", "histogram_pruning_equals_brute_force_40_comparisons", "7NXB_H3_span_regression", "heavy_fasta_parser", "failed_numbering_is_not_heavy", "engineered_domain_not_rejected_for_anchor_mutation"]})
    print("self-test passed", flush=True)


def analyze_simplefold():
    summaries = {}
    accessions = {}
    for name in ("swissprot_list.csv", "afesm_list.csv"):
        lines = [line.strip() for line in (SRC / name).read_text().splitlines() if line.strip()]
        assert all("," not in line and " " not in line for line in lines), "Expected headerless ID list"
        ids = {line.upper() for line in lines}
        accs = {m.group(1) for line in ids if (m := re.fullmatch(r"AF-(.+)-F\d+-MODEL_V\d+", line))}
        accessions[name] = accs
        summaries[name] = {"rows": len(lines), "unique_ids": len(ids), "afdb_accessions": len(accs),
            "non_afdb_ids": len(ids) - len(accs), "sha256": digest(SRC / name),
            "sequence_column_present": False}
    extended = json.loads(gzip.open(SRC / "afesme_dict.json.gz", "rt").read())
    assert all(isinstance(v["members"], list) for v in extended.values())
    summaries["afesme_dict.json"] = {"clusters": len(extended), "member_occurrences": sum(len(v["members"]) for v in extended.values()),
        "sha256": hashlib.sha256(gzip.open(SRC / "afesme_dict.json.gz", "rb").read()).hexdigest(), "role": "3B extended set; not substituted for the published 100M AFESM representatives"}
    del extended
    tree = json.loads((SRC / "simplefold_tree.json").read_text())
    public_paths = [r["path"] for r in tree["tree"]]
    save_json(OUT / "simplefold_repository_paths.json", public_paths)
    benchmark = json.loads((OUT / "benchmark_rcsb_entities.json").read_text())
    targets = json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())
    rows = []
    for t in targets:
        # SAbDab chain labels need not equal current author IDs; require sequence correspondence.
        es = [e for e in benchmark["entities"] if e["pdb"] == t["pdb"].upper()
              and t["hseq"] in "".join(e["entity_poly"]["pdbx_seq_one_letter_code_can"].split())]
        refs = sorted({r["database_accession"].upper() for e in es
            for r in e["rcsb_polymer_entity_container_identifiers"].get("reference_sequence_identifiers") or []
            if r["database_name"] == "UniProt"})
        hit = {name: sorted(set(refs) & accs) for name, accs in accessions.items()}
        rows.append({"tid": t["pdb"] + t["H"], "mapped_polymer_entity_ids": ";".join(e["rcsb_id"] for e in es),
            "uniprot_accessions": ";".join(refs), "swissprot_list_accession_hits": ";".join(hit["swissprot_list.csv"]),
            "afesm_list_accession_hits": ";".join(hit["afesm_list.csv"]),
            "interpretation": "exact frozen VH substring in entity; reference-accession screen only; missing mapping or no match does not exclude sequence homology"})
    save_csv(OUT / "simplefold_benchmark_identifier_screen.csv", rows)
    text = (SRC / "simplefold_paper.html").read_text()
    paragraphs = [html.unescape(re.sub("<[^>]+>", " ", s)) for s in re.findall(r"<p[^>]*>(.*?)</p>", text, re.S)]
    selected = [s for s in paragraphs if "SwissProt" in s and ("All SimpleFold models except" in s or "All models are trained with full" in s)]
    assert len(selected) >= 1
    save_json(OUT / "simplefold_recipe_evidence.json", {"source": "https://arxiv.org/html/2509.18480v1", "paragraphs": selected,
        "scope": "Published model-family recipe; exact historical benchmark checkpoint membership remains unverified"})
    summary = {"manifests": summaries, "github_tree_sha": tree["sha"], "repository_tree_truncated": tree["truncated"],
        "published_experimental_manifest_found_in_checked_tree": any("manifest" in p.lower() and p.endswith((".json", ".csv", ".txt")) for p in public_paths),
        "benchmark_targets": len(rows), "targets_with_chain_entity_mapping": sum(bool(r["mapped_polymer_entity_ids"]) for r in rows),
        "targets_with_uniprot_reference": sum(bool(r["uniprot_accessions"]) for r in rows),
        "targets_with_manifest_accession_match": sum(bool(r["swissprot_list_accession_hits"] or r["afesm_list_accession_hits"]) for r in rows),
        "limitations": ["No exact full experimental PDB training membership manifest found in the checked official README/tree; configuration references a local manifest.json, not a downloadable membership list.",
            "The recovered distilled manifests contain identifiers, not sequences. The accession screen is not a global sequence-similarity screen of approximately two million training structures.",
            "Absence of a benchmark UniProt reference is missing coverage, not a negative match.",
            "ESM2 sequence pretraining and exact historical checkpoint-specific memberships remain unaudited."]}
    save_json(OUT / "simplefold_manifest_summary.json", summary)
    print(json.dumps(summary), flush=True)


def cohort_summary(ids, scores, meta):
    import numpy as np
    models = ("SimpleFold", "IgFold", "NanoBodyBuilder2")
    result = {"n": len(ids), "n_pdbs": len({meta[t]["pdb"] for t in ids}), "tids": ids,
              "medians": {}, "paired_effects": {}, "rank_by_median": {}}
    if not ids:
        return result
    groups = {}
    for i, tid in enumerate(ids):
        groups.setdefault(meta[tid]["pdb"], []).append(i)
    blocks = list(groups.values())
    rng = np.random.default_rng(20261016)
    draws = rng.integers(len(ids), size=(4000, len(ids)))
    block_rng = np.random.default_rng(20261016)
    block_draws = [np.concatenate([blocks[j] for j in block_rng.integers(len(blocks), size=len(blocks))]) for _ in range(4000)]
    for model, s in zip(models, scores):
        result["medians"][model] = {key: float(np.median([s[t][key] for t in ids]))
            for key in ("chain_cdr1", "chain_cdr2", "chain_cdr3", "chain_framework", "fw_cdr1", "fw_cdr2", "fw_cdr3", "fw_framework")}
    for mode in ("chain", "fw"):
        result["rank_by_median"][mode] = sorted(models, key=lambda m: result["medians"][m][f"{mode}_cdr3"])
        for j, model in ((1, "IgFold"), (2, "NanoBodyBuilder2")):
            v = np.array([scores[0][t][f"{mode}_cdr3"] - scores[j][t][f"{mode}_cdr3"] for t in ids])
            result["paired_effects"][f"SF minus {model}: {mode}"] = {"median": float(np.median(v)),
                "ci95": np.quantile(np.median(v[draws], axis=1), [.025, .975]).tolist(),
                "pdb_block_ci95": np.quantile([np.median(v[k]) for k in block_draws], [.025, .975]).tolist()}
    return result


def validate_manifest_accessions():
    targets = {t["pdb"] + t["H"]: t for t in json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())}
    flagged = [r for r in read_csv(OUT / "simplefold_benchmark_identifier_screen.csv")
               if r["swissprot_list_accession_hits"] or r["afesm_list_accession_hits"]]
    aligner = make_aligner()
    rows = []
    for row in flagged:
        accessions = sorted(set((row["swissprot_list_accession_hits"] + ";" + row["afesm_list_accession_hits"]).strip(";").split(";")))
        for accession in accessions:
            path = fetch(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta", f"uniprot_{accession}.fasta")
            lines = path.read_text().splitlines()
            assert lines[0].startswith(">")
            canonical = "".join(lines[1:])
            vh = targets[row["tid"]]["hseq"]
            rows.append({"tid": row["tid"], "accession": accession, "source_header": lines[0],
                "canonical_length": len(canonical), "benchmark_vh_length": len(vh),
                "exact_vh_substring_of_canonical": vh in canonical,
                "global_vh_identity_to_full_canonical": identity(aligner, vh, canonical),
                "interpretation": "current canonical sequence check, not historical AFDB-model sequence verification"})
    save_csv(OUT / "simplefold_accession_sequence_checks.csv", rows)
    print(json.dumps(rows, indent=2), flush=True)


def summarize():
    import numpy as np
    old = json.loads((OLD / "training_overlap_audit.json").read_text())
    targets = json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())
    meta = {t["pdb"] + t["H"]: t for t in targets}
    baseline = [{r["tid"]: r for r in read_csv(OLD / f"{name}_nearest_training.csv")} for name in
                ("nbb2", "igfold", "simplefold_pre_cutoff_proxy")]
    extensions = {name: {r["tid"]: r for r in read_csv(OUT / f"{name}_nearest.csv")}
        for name in ("nbb2_extension", "af_oas_paired", "af_oas_unpaired") if (OUT / f"{name}_nearest.csv").exists()}
    assert "nbb2_extension" in extensions
    scores = [{r["tid"]: r for r in json.loads((ROOT / p).read_text())} for p in (
        "outputs/structure/simplefold/metrics/per_target_scores.json",
        "outputs/structure/igfold/metrics/per_target_scores_igfold.json",
        "outputs/structure/nbb2/per_target_scores_nbb2.json")]
    paired = sorted(set.intersection(*(set(s) for s in scores)))
    assert len(paired) == 588
    combined, changed = [], []
    coverage_counts = {}
    for name, block in extensions.items():
        coverage_counts[name] = {"all597": {}, "paired588": {}}
        for label, ids in (("all597", meta), ("paired588", paired)):
            coverage_counts[name][label] = {"vhh_ge90": sum(float(block[t]["max_vhh_identity"]) >= .9 for t in ids),
                "cdr3_ge80": sum(float(block[t]["max_cdr3_identity"]) >= .8 for t in ids)}
    for tid in paired:
        old_vh = max(float(b[tid]["max_vhh_identity"]) for b in baseline)
        old_h3 = max(float(b[tid]["max_cdr3_identity"]) for b in baseline)
        vh = max([old_vh] + [float(b[tid]["max_vhh_identity"]) for b in extensions.values()])
        h3 = max([old_h3] + [float(b[tid]["max_cdr3_identity"]) for b in extensions.values()])
        no_pdb = all(b[tid]["pdb_overlap"] == "False" for b in baseline)
        common_atoms = len({s[tid]["n_atoms_matched"] for s in scores}) == 1
        old_pass = no_pdb and old_vh < .9 and old_h3 < .8 and tid != "9RITB" and common_atoms
        new_pass = no_pdb and vh < .9 and h3 < .8 and tid != "9RITB" and common_atoms
        r = {"tid": tid, "no_known_training_pdb": no_pdb, "old_max_vhh_identity": old_vh,
            "old_max_cdr3_identity": old_h3, "extended_max_vhh_identity": vh, "extended_max_cdr3_identity": h3,
            "same_matched_atom_count": common_atoms, "old_pass_before_clustering": old_pass,
            "extended_pass_before_clustering": new_pass,
            "extension_vhh_threshold_sources": ";".join(n for n, b in extensions.items() if float(b[tid]["max_vhh_identity"]) >= .9),
            "extension_cdr3_threshold_sources": ";".join(n for n, b in extensions.items() if float(b[tid]["max_cdr3_identity"]) >= .8)}
        combined.append(r)
        if old_pass and not new_pass:
            changed.append(r)
    save_csv(OUT / "extended_filter_flags.csv", combined)
    if changed:
        save_csv(OUT / "newly_excluded_targets.csv", changed)
    cache = np.load(OLD / "pairwise_identity.npz")
    positions = {t: i for i, t in enumerate(cache["tids"].tolist())}
    vhh_matrix = cache["vh"]
    def representatives(ids):
        parent = {t: t for t in ids}
        def find(t):
            while parent[t] != t:
                parent[t] = parent[parent[t]]
                t = parent[t]
            return t
        for i, tid in enumerate(ids):
            for other in ids[:i]:
                if vhh_matrix[positions[tid], positions[other]] >= .9:
                    parent[find(tid)] = find(other)
        groups = {}
        for t in ids:
            groups.setdefault(find(t), []).append(t)
        return sorted(min(g, key=lambda t: (float(meta[t]["resolution"]), t)) for g in groups.values())
    old_pass = [r["tid"] for r in combined if r["old_pass_before_clustering"]]
    assert len(old_pass) == 403
    old_ids = representatives(old_pass)
    assert old_ids == old["common_lower_risk_subset"]["tids"] and len(old_ids) == 398
    new_pass = [r["tid"] for r in combined if r["extended_pass_before_clustering"]]
    new_ids = representatives(new_pass)
    target_rows = [{"tid": tid, "pdb": meta[tid]["pdb"], "chain": meta[tid]["H"],
                    "instance": meta[tid]["instance"], "resolution": meta[tid]["resolution"]} for tid in new_ids]
    save_csv(OUT / "extended_lower_risk_targets.csv", target_rows)
    pdb_excluded = [r["tid"] for r in combined if r["no_known_training_pdb"]]
    assert len(pdb_excluded) == 453
    complete = all(n in extensions for n in ("nbb2_extension", "af_oas_paired", "af_oas_unpaired"))
    report = {"status": "FEASIBLE_SOURCE_EXTENSION_COMPLETE_B2_OVERALL_PARTIAL" if complete else "INTERIM_NBB2_COMPLETE_OAS_PENDING",
        "not_certified_training_independent": True, "extension_pools": list(extensions),
        "extension_flag_counts": coverage_counts, "before_clustering_old": len(old_pass),
        "before_clustering_extended": len(new_pass), "newly_excluded_before_clustering": [r["tid"] for r in changed],
        "removed_old_representatives": sorted(set(old_ids) - set(new_ids)),
        "added_representatives": sorted(set(new_ids) - set(old_ids)), "cohorts": {}}
    for label, ids in (("original_common", paired), ("known_training_pdb_excluded", pdb_excluded),
                       ("16oct_lower_risk", old_ids), ("extended_lower_risk", new_ids)):
        report["cohorts"][label] = cohort_summary(ids, scores, meta)
    for model in ("SimpleFold", "IgFold", "NanoBodyBuilder2"):
        for mode in ("chain_cdr3", "fw_cdr3"):
            assert report["cohorts"]["16oct_lower_risk"]["medians"][model][mode] == old["lower_risk_nonredundant"]["medians"][model][mode]
    save_json(OUT / "overlap_extension_summary.json", report)
    if not complete:
        save_json(OUT / "nbb2_only_cohort_summary.json", report)
        save_csv(OUT / "nbb2_only_lower_risk_targets.csv", target_rows)
    save_csv(OUT / "extended_lower_risk_scores.csv", [{"tid": tid, "pdb": meta[tid]["pdb"], "model": model,
        **{key: s[tid][key] for key in ("chain_cdr1", "chain_cdr2", "chain_cdr3", "chain_framework", "fw_cdr1", "fw_cdr2", "fw_cdr3", "fw_framework")}}
        for tid in new_ids for model, s in zip(("SimpleFold", "IgFold", "NanoBodyBuilder2"), scores)])
    if complete:
        nbb = json.loads((OUT / "nbb2_recovery_summary.json").read_text())
        oas = {name: json.loads((OUT / f"{name}_coverage.json").read_text())
               for name in ("af_oas_paired", "af_oas_unpaired")}
        lines = ["# B2 source-bound overlap extension — 18 October 2026", "",
            "Status: feasible public-source extension complete; full model-specific exposure audit remains PARTIAL. This is not a training-independent benchmark certification.", "",
            "## Scope and validation", "",
            f"Recovered heavy-domain candidates for all {nbb['pdbs_entity_sequences_recovered']} previously unmapped NanoBodyBuilder2 manifest PDBs: {nbb['accepted_heavy_domains']} accepted domains and {nbb['unique_vh_h3_pairs']} distinct VH/H3 pairs. ANARCII {nbb['anarcii_version']} numbered every recovered protein entity rather than assuming a chain from its name. Benchmark H3 validation agreed on 597/597; full-domain validation agreed on 593/597 because of terminal trimming. Candidate recovery does not prove historical training-chain membership.", "",
            "IgFold source: https://doi.org/10.5281/zenodo.7820263. Both complete OAS ZIPs passed their official published MD5 checksums at collection. Coverage and archive SHA256 values are in the per-archive coverage JSON and source receipts.", ""]
        for name, block in oas.items():
            lines += [f"- {name}: {block['fasta_files']:,} legitimate FASTAs, comprising {block['heavy_pairs']:,} heavy records and {block['light_only_fasta_members']:,} light-only records excluded from VHH screening; {block['unique_heavy_h3_pairs']:,} distinct VH/H3 pairs; {block['ignored_resource_fork_fasta_members']} Apple resource-fork FASTA members excluded; {len(block['unavailable_h3_members'])} unresolved H3 records. Coordinate-member status: {block['coordinate_status_counts']}. H3 extraction: {block['h3_method_counts']}."]
        lines += ["", "HTML error pages packaged as PDB files were not treated as coordinates. Their valid heavy FASTA sequences were retained; H3 recovery used ANARCII 2.0.8 IMGT 105–117 and is explicitly labelled. Unknown H3s, if any, are missing coverage, not negative matches.", "",
            "The nearest-neighbour screen uses the frozen global affine-gap definition (+2/-1/-5/-0.5; identities divided by all alignment columns). An amino-acid histogram bound safely prunes candidates without changing the maximum. The self-test compares pruned results with 40 brute-force comparisons and checks the saved 7NXB H3-span regression. Separate VH/H3 threshold flags are not added as independent counts.", "",
            "## Same-score sensitivity comparison", "", "Model order: SimpleFold, IgFold, NanoBodyBuilder2. Medians below use retained prediction scores; no folding model was run.", ""]
        for label, block in report["cohorts"].items():
            medians = block["medians"]
            lines += [f"- {label}: n={block['n']}, PDB entries={block['n_pdbs']}; chain H3 medians " + ", ".join(f"{medians[m]['chain_cdr3']:.4f}" for m in ("SimpleFold", "IgFold", "NanoBodyBuilder2")) + "; framework H3 medians " + ", ".join(f"{medians[m]['fw_cdr3']:.4f}" for m in ("SimpleFold", "IgFold", "NanoBodyBuilder2")) + " Å."]
        lines += ["", f"Candidates before clustering: initial {len(old_pass)}, expanded {len(new_pass)}. Representatives: initial 398, expanded {len(new_ids)}. Newly failing candidate IDs: " + ", ".join(report["newly_excluded_before_clustering"]) + ". Resolution/identifier-based representative selection never uses errors. Paired target-bootstrap and PDB-block intervals (4,000 draws, seed 20261016) are retained in overlap_extension_summary.json; thresholds are exploratory, not preregistered.", "",
            "Filtering changes target composition and exposure risk together. Ranking differences do not establish leakage as their sole cause or universal superiority of any model.", "",
            "## SimpleFold and unresolved coverage", "",
            "The official SwissProt and AFESM files contain 269,003 and 2,036,414 IDs, not sequence columns. The extended dictionary has 8,616,595 member occurrences and belongs to the extended recipe, not the 100M representative list. Checked official GitHub README/tree/config sources do not supply the full experimental membership manifest. Source receipts, recipe excerpts and the repository tree SHA are retained.", "",
            "Current RCSB entity mapping covers 594/597 VHH segments. Six entity-level manifest-accession flags were checked against current UniProt FASTAs: none contains the corresponding frozen VHH substring, and the accession descriptions are bacterial proteins. These flags are not accepted as confirmed VHH overlap. This current canonical check is not historical AFDB-model sequence verification.", "",
            "Unresolved: full SimpleFold experimental/distilled sequence-level membership/similarity, language-model sequence pretraining, exact historical checkpoint exposure, and antigen/epitope/interface similarity. Available published OAS sources have now been screened, but do not prove the exact historical IgFold training snapshot. Missing evidence is disclosed, not interpreted as independence.", "",
            "## Reproduction", "",
            "Use isolated uv with Python 3.12, Biopython 1.87, NumPy 2.2.6 and ANARCII 2.0.8. Stages: collect-rcsb, collect-obsolete, collect-benchmark, collect-oas, collect-simplefold, number-nbb, compare-nbb, compare-oas, analyze-simplefold, validate-simplefold-accessions, self-test, summarize. Existing validated source caches are reused. Collection requires network access; full pairwise screening requires CPU time, not a folding GPU. The frozen_input_sha256.json guard checks all original evidence files before and after each stage.", ""]
        (OUT / "REPORT.md").write_text("\n".join(lines))
    print(json.dumps({"status": report["status"], "n": len(new_ids), "new_exclusions": len(changed),
        "medians": report["cohorts"]["extended_lower_risk"]["medians"], "ranks": report["cohorts"]["extended_lower_risk"]["rank_by_median"]}), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["collect-rcsb", "collect-oas", "collect-simplefold", "collect-obsolete", "collect-benchmark", "number-nbb", "compare-nbb", "compare-oas", "self-test", "analyze-simplefold", "validate-simplefold-accessions", "summarize"])
    args = p.parse_args()
    SRC.mkdir(parents=True, exist_ok=True)
    frozen_inputs()
    {"collect-rcsb": collect_rcsb, "collect-oas": collect_oas,
     "collect-simplefold": collect_simplefold, "collect-obsolete": collect_obsolete,
     "collect-benchmark": lambda: collect_rcsb(benchmark=True), "number-nbb": number_nbb,
     "compare-nbb": compare_nbb, "compare-oas": compare_oas, "self-test": self_test,
     "analyze-simplefold": analyze_simplefold, "validate-simplefold-accessions": validate_manifest_accessions,
     "summarize": summarize}[args.stage]()
    frozen_inputs()


if __name__ == "__main__":
    main()
