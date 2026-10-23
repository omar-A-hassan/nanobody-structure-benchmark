"""Fetch public source metadata, not model predictions. No journal or archive uploads."""
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
import csv
import json
import tarfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/revision_2026-10-16/sources"
OUT.mkdir(parents=True, exist_ok=True)

# Stop the archive stream after the metadata; do not download all coordinate files.
names = {"README.md", "ab_split_sd.csv", "ab_split.csv"}
if not all((OUT / n).exists() for n in names):
    url = "https://zenodo.org/records/20083995/files/splits.tar.gz?download=1"
    with urlopen(url, timeout=90) as response:
        with tarfile.open(fileobj=response, mode="r|gz") as archive:
            for member in archive:
                leaf = Path(member.name).name
                if leaf in names and member.isfile():
                    payload = archive.extractfile(member).read()
                    assert payload
                    (OUT / leaf).write_bytes(payload)
                    print("Fetched", member.name, len(payload), flush=True)
                if all((OUT / n).exists() for n in names):
                    break
assert all((OUT / n).exists() for n in names), "Source archive lacks expected metadata"

# Initial public release is not the same as deposition date saved in targets_vhh.json.
targets = json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())
ids = sorted({t["pdb"].upper() for t in targets})
outfile = OUT / "rcsb_release_dates.json"
found = json.loads(outfile.read_text()) if outfile.exists() else {}
for start in range(0, len(ids), 50):
    batch = [p for p in ids[start:start + 50] if p not in found]
    if not batch:
        continue
    query = "{entries(entry_ids:" + json.dumps(batch) + "){rcsb_id rcsb_accession_info{initial_release_date deposit_date}}}"
    url = "https://data.rcsb.org/graphql?" + urlencode({"query": query})
    with urlopen(url, timeout=90) as response:
        data = json.load(response)
    assert not data.get("errors"), data.get("errors")
    for entry in data["data"]["entries"]:
        assert entry is not None
        found[entry["rcsb_id"].upper()] = entry["rcsb_accession_info"]
    outfile.write_text(json.dumps(found, indent=2) + "\n")
    print("RCSB release dates", len(found), "/", len(ids), flush=True)
missing = set(ids) - set(found)
if missing == {"9RIT"}:
    url = "https://www.rcsb.org/structure/removed/9RIT"
    with urlopen(url, timeout=90) as response:
        html = response.read().decode()
    assert all(s in html for s in ("2025-11-05", "2025-06-11", "2026-09-02", "9S9E"))
    found["9RIT"] = {"initial_release_date": "2025-11-05", "deposit_date": "2025-06-11",
                     "status": "obsolete", "obsolete_date": "2026-09-02", "replacement": "9S9E", "source": url}
    outfile.write_text(json.dumps(found, indent=2) + "\n")
assert set(ids) <= set(found)
with (OUT.parent / "benchmark_release_dates.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["tid", "pdb", "source_deposition_date", "public_release_date",
                                           "pdb_status", "replacement", "released_after_simplefold_cutoff"])
    writer.writeheader()
    for t in targets:
        info = found[t["pdb"].upper()]
        released = info["initial_release_date"][:10]
        writer.writerow({"tid": t["pdb"] + t["H"], "pdb": t["pdb"], "source_deposition_date": t["PDBdepo"],
                         "public_release_date": released, "pdb_status": info.get("status", "current"),
                         "replacement": info.get("replacement", ""), "released_after_simplefold_cutoff": released > "2020-05-01"})
print("Verified all release dates; 9RIT is obsolete and must be flagged", flush=True)
