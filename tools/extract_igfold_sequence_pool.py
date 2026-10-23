"""Extract the exact audited IgFold experimental VH/H3 pool from its verified source ZIP."""
from pathlib import Path
import argparse, hashlib, json, zipfile
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("archive",type=Path)
p.add_argument("output",type=Path)
a=p.parse_args()
zpath = a.archive
assert hashlib.md5(zpath.read_bytes()).hexdigest() == "5a13e097ac6dfd03c8f3fcf93f6eb065"
igseq, igpdbs = set(), set()
with zipfile.ZipFile(zpath) as z:
    for name in z.namelist():
        if not name.endswith(".fasta"):
            continue
        pdb = Path(name).stem.upper()
        igpdbs.add(pdb)
        heavy = ""
        active = False
        for line in z.read(name).decode().splitlines():
            if line.startswith(">"):
                active = line.rsplit(":", 1)[-1] == "H"
            elif active:
                heavy += line.strip()
        if not heavy:
            continue
        # Chothia H93–H102 corresponds to the IMGT H105–H117 span used in this benchmark.
        seen, h3 = set(), ""
        for line in z.read(str(Path(name).with_suffix(".pdb"))).decode().splitlines():
            if line.startswith("ATOM") and line[21] == "H" and line[12:16].strip() == "CA":
                key = line[22:27]
                if key not in seen and 93 <= int(line[22:26]) <= 102:
                    aa3 = line[17:20]
                    aa = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
                          "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
                          "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}[aa3]
                    h3 += aa
                    seen.add(key)
        igseq.add((heavy, h3))
        if pdb == "7NXB":
            assert h3 == "AKDGGKLWVYYFDY", "IMGT-equivalent H3 span must include Chothia H93/H94"
assert len(igpdbs) == 4275

a.output.write_text(json.dumps({"heavy_h3_pairs":sorted(igseq),"pdb_ids":sorted(igpdbs)},indent=2)+"\n")
