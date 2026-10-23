"""Assemble a public analysis snapshot from an explicit local source selection.

Run in an isolated uv environment with NumPy. No network or inference is used.
Existing destination artifacts must match; this never cleans a source directory.
"""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import zipfile

import numpy as np

DEST = Path(__file__).resolve().parents[1]
SECRET = [r"hf_[A-Za-z0-9]{20,}", r"KGAT_[A-Za-z0-9_]{12,}",
          r"gh[pousr]_[A-Za-z0-9]{20,}", r"wandb_[A-Za-z0-9]{20,}"]
TEXT = {".py", ".json", ".ipynb", ".csv", ".md", ".txt", ".log", ".yaml", ".fasta", ".html"}


def digest(path, algorithm="sha256"):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def write_regression_fixture(source):
    archive_path = source / "outputs/revision_2026-10-16/sources/igfold_sabdab.zip"
    assert digest(archive_path, "md5") == "5a13e097ac6dfd03c8f3fcf93f6eb065"
    with zipfile.ZipFile(archive_path) as archive:
        member = next(n for n in archive.namelist() if n.lower().endswith("7nxb.pdb"))
        payload = archive.read(member)
    target = DEST / "data/regression/igfold_7nxb.pdb"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(payload)
    save(target.with_suffix(".receipt.json"), {"source_archive_md5": digest(archive_path, "md5"),
        "source_archive_sha256": digest(archive_path), "archive_member": member, "sha256": digest(target),
        "source_url": "https://zenodo.org/records/7820263/files/sabdab.zip?download=1", "transformation": "Unmodified archive member"})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--figure-source", type=Path, required=True)
    parser.add_argument("--fixtures-only", action="store_true")
    parser.add_argument("--refresh-receipts", action="store_true")
    parser.add_argument("--sync-latest-code", action="store_true")
    args = parser.parse_args()
    source = args.source.resolve()
    assert source != DEST and not DEST.is_relative_to(source)
    if args.sync_latest_code:
        assert not (DEST / "inference/repeat_outputs").exists(), "Coordinate package updates after repeat results are added"
        selection_path = DEST / "provenance/source_selection.json"
        selection = json.loads(selection_path.read_text())
        names = ["src/make_revision_workflow.py", "src/verify_archive_reruns.py"]
        names += [f"kaggle/kernels/revision-archive-20261023/{m}/kernel-metadata.json" for m in ["simplefold", "nbb2"]]
        for name in names:
            source_name = name
            if name.endswith("kernel-metadata.json"):
                model = Path(name).parent.name
                source_name = f"outputs/revision_2026-10-23/repeat_runs/{model}/source/kernel-metadata.json"
            path = source / source_name
            target = DEST / name
            assert path.is_file()
            note = "byte-identical latest source; notebooks were not regenerated"
            if name.endswith("kernel-metadata.json"):
                metadata = json.loads(path.read_text())
                metadata["code_file"] = "repeat.ipynb"
                save(target, metadata)
                note = "Latest locally recovered actual remote metadata; code_file adjusted to repeat.ipynb; notebook not regenerated"
            else:
                shutil.copyfile(path, target)
            records = [r for r in selection["files"] if r["path"] == name]
            row = records[0] if records else {}
            row.update(path=name, source=source_name, source_root="study", source_sha256=digest(path),
                       packaged_sha256=digest(target), transformation=note)
            if not records:
                selection["files"].append(row)
        save(selection_path, selection)
        return
    if args.refresh_receipts:
        selection_path = DEST / "provenance/source_selection.json"
        selection = json.loads(selection_path.read_text())
        for name in ["source_manifest.json", "output_manifest.json"]:
            old = f"outputs/revision_2026-10-17/reference_reconstruction/{name}"
            new = f"provenance/reference_reconstruction_historical_{name}"
            if (DEST / old).exists():
                assert not (DEST / new).exists()
                (DEST / old).rename(DEST / new)
            for row in selection["files"]:
                if row["path"] == old:
                    row["path"] = new
                    row["transformation"] = "Byte-identical historical receipt, relocated to distinguish its scope from the public inventory"
        for row in selection["files"]:
            path = DEST / row["path"]
            assert path.is_file(), row["path"]
            row["packaged_sha256"] = digest(path)
            if row["source_sha256"] != row["packaged_sha256"] and row["transformation"] == "byte-identical":
                row["transformation"] = "Public storage, path or report-neutrality adaptation; see docs/DEPENDENCIES.md"
        selection["derived_fixtures"] = [{"path": "data/regression/igfold_7nxb.pdb", "receipt": "data/regression/igfold_7nxb.receipt.json"}]
        save(selection_path, selection)
        return
    if args.fixtures_only:
        write_regression_fixture(source)
        return
    records = []

    def copy(relative, destination=None, root=source):
        path = root / relative
        target = DEST / (destination or relative)
        assert path.is_file() and not path.is_symlink(), path
        assert target.is_relative_to(DEST) and ".." not in target.parts
        if path.suffix in TEXT:
            assert not any(re.search(p, path.read_text()) for p in SECRET), relative
        target.parent.mkdir(parents=True, exist_ok=True)
        note = "byte-identical"
        if path.suffix == ".ipynb":
            nb = json.loads(path.read_text())
            for cell in nb["cells"]:
                if cell["cell_type"] == "code":
                    cell["outputs"] = []
                    cell["execution_count"] = None
                cell.pop("attachments", None)
            nb["metadata"].pop("widgets", None)
            target.write_text(json.dumps(nb, indent=1) + "\n")
            note = "notebook outputs/execution counts/widgets removed; source cells preserved"
        elif target.exists():
            assert digest(path) == digest(target), f"Existing artifact differs: {target}"
        else:
            shutil.copyfile(path, target)
        records.append({"path": str(target.relative_to(DEST)), "source": str(relative),
                        "source_root": "study" if root == source else "original_figure_source",
                        "source_sha256": digest(path), "packaged_sha256": digest(target),
                        "transformation": note})
        return target

    def compressed(relative):
        path = source / relative
        target = DEST / (relative + ".gz")
        target.parent.mkdir(parents=True, exist_ok=True)
        with path.open("rb") as inp, target.open("wb") as out:
            with gzip.GzipFile(filename="", mode="wb", fileobj=out, mtime=0) as stream:
                shutil.copyfileobj(inp, stream)
        with gzip.open(target, "rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == digest(path)
        records.append({"path": str(target.relative_to(DEST)), "source": relative,
                        "source_root": "study", "source_sha256": digest(path),
                        "packaged_sha256": digest(target), "transformation": "lossless gzip, mtime=0"})

    scripts = ["reviewer_revision_analysis", "revision_source_metadata", "revision_training_audit",
               "revision_overlap_extension", "reconstruct_revision_references", "revision_reference_sensitivity",
               "verify_retained_structure_scoring", "esm_probe_analyze", "reproduce_esm_probe",
               "recreate_reused_figures", "make_revision_workflow", "verify_archive_reruns", "replay_revision_nbb2"]
    for name in scripts:
        copy(f"src/{name}.py")
    copy("src/consolidate_stats.py", "provenance/original_scripts/consolidate_stats.py")
    for path in sorted((source / "src/tnp").glob("*.py")):
        copy(str(path.relative_to(source)))
    for path in sorted((source / "data/inputs").glob("*.json")):
        copy(str(path.relative_to(source)))
    for folder in ["structure", "hard", "hitrate", "trajectory", "tnp", "ensemble", "visualization"]:
        for path in sorted((source / "outputs" / folder).rglob("*")):
            if path.is_file() and path.suffix in {".json", ".csv", ".npy", ".pdb", ".log"}:
                copy(str(path.relative_to(source)))
    for path in sorted((source / "outputs/esm-probe").glob("*")):
        if path.suffix in {".json", ".log"}:
            copy(str(path.relative_to(source)))

    old = source / "outputs/revision_2026-10-16"
    for path in sorted(old.iterdir()):
        if path.is_file() and path.suffix in {".csv", ".json", ".npz", ".fasta"} and path.name != "integrity_manifest.json":
            copy(str(path.relative_to(source)))
        elif path.suffix in {".png", ".pdf"}:
            copy(str(path.relative_to(source)), "figures/saved/" + path.name)
    for name in ["README.md", "ab_split_sd.csv", "nbb2_nanobody_data.csv", "rcsb_release_dates.json"]:
        copy("outputs/revision_2026-10-16/sources/" + name)
    compressed("outputs/revision_2026-10-16/sources/ab_split.csv")

    # Derive exactly the audited experimental pool without redistributing a 149 MB ZIP.
    audit = (source / "src/revision_training_audit.py").read_text()
    block = audit[audit.index('zpath = SRC / "igfold_sabdab.zip"'):audit.index('\naligner = PairwiseAligner')]
    scope = {"SRC": old / "sources", "hashlib": hashlib, "zipfile": zipfile, "Path": Path}
    exec(compile(block, "verified_igfold_pool_extraction", "exec"), scope)
    assert len(scope["igseq"]) == 2649 and len(scope["igpdbs"]) == 4275
    pool = {"heavy_h3_pairs": sorted(scope["igseq"]), "pdb_ids": sorted(scope["igpdbs"]),
            "source_archive": {"url": "https://zenodo.org/records/7820263/files/sabdab.zip?download=1",
                               "md5": digest(scope["zpath"], "md5"), "sha256": digest(scope["zpath"]),
                               "bytes": scope["zpath"].stat().st_size},
            "extraction": "Original audit block; Chothia H93–H102 including insertion codes; 7NXB regression asserted"}
    save(DEST / "outputs/revision_2026-10-16/sources/igfold_experimental_sequence_pool.json", pool)
    copy("src/revision_training_audit.py", "provenance/original_scripts/revision_training_audit.py")
    extractor = '"""Extract the exact audited IgFold experimental VH/H3 pool from its verified source ZIP."""\n'
    extractor += "from pathlib import Path\nimport argparse, hashlib, json, zipfile\n"
    extractor += 'p=argparse.ArgumentParser(description=__doc__)\np.add_argument("archive",type=Path)\np.add_argument("output",type=Path)\na=p.parse_args()\n'
    extractor += block.replace('zpath = SRC / "igfold_sabdab.zip"', 'zpath = a.archive')
    extractor += '\na.output.write_text(json.dumps({"heavy_h3_pairs":sorted(igseq),"pdb_ids":sorted(igpdbs)},indent=2)+"\\n")\n'
    (DEST / "tools/extract_igfold_sequence_pool.py").write_text(extractor)

    for folder in ["overlap_extension", "reference_reconstruction", "reference_sensitivity", "reused_figures"]:
        base = source / "outputs/revision_2026-10-17" / folder
        for path in sorted(base.rglob("*")):
            if not path.is_file() or path.name.startswith("partial_"):
                continue
            if path.suffix in {".log", ".part", ".bin"} or path.name in {"REPORT.md", "frozen_input_sha256.json"}:
                continue
            if folder == "reference_reconstruction" and "sources" in path.relative_to(base).parts:
                if path.suffix != ".cif" and path.name not in {"download_receipt.json"}:
                    continue
            if path.name in {"af_oas_paired.zip", "af_oas_unpaired.zip", "splits.tar.gz"}:
                continue
            if path.name == "afesme_dict.json":
                compressed(str(path.relative_to(source)))
            elif path.stat().st_size < 100_000_000:
                copy(str(path.relative_to(source)))
    # Small, separately labelled inference check; weights and native executables are excluded.
    base = source / "outputs/revision_2026-10-17/nbb2_regeneration_check"
    for path in sorted(base.rglob("*")):
        if path.is_file() and path.name != "failed_packaged_hmmer_wrapper.json" and path.suffix in {".json", ".csv", ".py", ".pdb", ".lock", ".txt"}:
            if "weights" not in path.relative_to(base).parts and "bin" not in path.relative_to(base).parts:
                copy(str(path.relative_to(source)))

    kernels = ["nbb2-bench", "nbb2-bench-v5", "sf-hard", "sf-hitrate", "sf-traj", "sf-traj2",
               "sf-traj-tfm-597", "sf-tfm-only-597", "sf-ensemble", "tnp-sf", "tnp-specialists",
               "esm-probe", "sf-viz", "nbb2-viz", "revision-archive-20261023"]
    for folder in kernels:
        for path in sorted((source / "kaggle/kernels" / folder).rglob("*")):
            if path.is_file() and path.suffix in {".json", ".ipynb"}:
                copy(str(path.relative_to(source)))
    for path in sorted((source / "kaggle/notebooks/legacy").glob("*.ipynb")):
        copy(str(path.relative_to(source)))
    for model in ["simplefold", "nbb2"]:
        relative = f"kaggle/kernels/revision-archive-20261023/{model}/kernel-metadata.json"
        meta_source = f"outputs/revision_2026-10-23/repeat_runs/{model}/source/kernel-metadata.json"
        meta = json.loads((source / meta_source).read_text())
        meta["code_file"] = "repeat.ipynb"
        save(DEST / relative, meta)
        next(r for r in records if r["path"] == relative)["transformation"] = "Actual locally recovered remote metadata; code_file adjusted to packaged repeat.ipynb"

    for path in sorted((source / "outputs/revision_2026-10-23/recovery").rglob("*")):
        if path.is_file() and path.suffix in {".json", ".ipynb", ".csv", ".log"}:
            copy(str(path.relative_to(source)))
    for path in sorted((source / "outputs/revision_2026-10-23/probe_reproduction").glob("*")):
        if path.suffix in {".json", ".npz"}:
            copy(str(path.relative_to(source)))
    copy(args.figure_source.name, "paper_figure_sources/make_figures.py", root=args.figure_source.parent)
    for stem in ["figure2_uncertainty", "figure3_prior_coverage", "figure4_esm_probe"]:
        for suffix in [".png", ".pdf"]:
            copy("figures/" + stem + suffix, "paper_figure_sources/original_assets/" + stem + suffix, root=args.figure_source.parent)
    handoff = "outputs/revision_2026-10-17/reproducibility_handoff/"
    for name in ["SimpleFold-LICENSE.txt", "ImmuneBuilder-LICENSE.txt", "sabdab2_record_metadata.json", "igfold_record_metadata.json", "reconstructed_worker_assignments.csv"]:
        copy(handoff + name, "provenance/" + name)
    history = json.loads((source / handoff / "historical_execution_provenance.json").read_text())
    history["igfold_inference"] = "Original inference package and execution provenance unavailable in the supplied artifacts"
    save(DEST / "provenance/historical_execution_provenance.json", history)

    embeddings = source / "outputs/esm-probe/esm_probe_embeddings.npz"
    with np.load(embeddings) as z:
        arrays = {k: z[k] for k in z.files}
    shards = []
    for start in range(0, 37, 8):
        end = min(start + 8, 37)
        target = DEST / f"data/esm_probe/embeddings_{start:02d}_{end:02d}.npz"
        target.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(target, **{k: v[:, start:end] for k, v in arrays.items()})
        with np.load(target) as z:
            assert all(np.array_equal(z[k], arrays[k][:, start:end]) for k in arrays)
        shards.append({"path": str(target.relative_to(DEST)), "start_layer": start, "end_layer_exclusive": end,
                       "sha256": digest(target)})
    save(DEST / "data/esm_probe/embedding_manifest.json", {
        "original_sha256": digest(embeddings), "original_bytes": embeddings.stat().st_size,
        "arrays": {k: {"shape": v.shape, "dtype": str(v.dtype), "sha256_c_order": hashlib.sha256(v.tobytes()).hexdigest()} for k, v in arrays.items()},
        "transformation": "Lossless layer slicing; no numeric conversion", "shards": shards})

    # Public adaptations: only storage, portable paths and the curated frozen receipt.
    target = DEST / "src/revision_training_audit.py"
    text = target.read_text().replace("import zipfile", "import zipfile\nimport gzip")
    text = text.replace('with (SRC / "ab_split.csv").open() as f:', 'with gzip.open(SRC / "ab_split.csv.gz", "rt", newline="") as f:')
    replacement = 'pool = json.loads((SRC / "igfold_experimental_sequence_pool.json").read_text())\nigseq = {tuple(row) for row in pool["heavy_h3_pairs"]}\nigpdbs = set(pool["pdb_ids"])\nassert len(igseq) == 2649 and len(igpdbs) == 4275\n'
    assert block in text
    text = text.replace(block, replacement)
    anchor = 'report["sources_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in\n        (SRC / "ab_split.csv", SRC / "nbb2_nanobody_data.csv", zpath)}'
    assert anchor in text
    text = text.replace(anchor, 'with gzip.open(SRC / "ab_split.csv.gz", "rb") as stream:\n    ab_hash = hashlib.file_digest(stream, "sha256").hexdigest()\nreport["sources_sha256"] = {"ab_split.csv": ab_hash, "nbb2_nanobody_data.csv": hashlib.sha256((SRC / "nbb2_nanobody_data.csv").read_bytes()).hexdigest(), "igfold_sabdab.zip": pool["source_archive"]["sha256"]}')
    target.write_text(text)
    target = DEST / "src/revision_overlap_extension.py"
    text = target.read_text()
    begin, end = text.index("def frozen_inputs():"), text.index("\n\ndef fetch(")
    text = text[:begin] + '''def frozen_inputs():
    """Verify the explicitly curated PUBLIC subset, including compressed sources."""
    receipt = OUT / "frozen_input_sha256.json"
    expected = json.loads(receipt.read_text())
    now = {name: digest(ROOT / name) for name in expected}
    assert now == expected, "Public frozen inputs changed; use a fresh checkout for replay"
    return now
''' + text[end:]
    text = text.replace('(SRC / "afesme_dict.json").read_text()', 'gzip.open(SRC / "afesme_dict.json.gz", "rt").read()')
    text = text.replace('digest(SRC / "afesme_dict.json")', 'hashlib.sha256(gzip.open(SRC / "afesme_dict.json.gz", "rb").read()).hexdigest()')
    old_test = '''    with zipfile.ZipFile(OLD / "sources/igfold_sabdab.zip") as archive:
        name = next(n for n in archive.namelist() if n.lower().endswith("7nxb.pdb"))
        assert pdb_h3(archive.read(name).decode()) == "AKDGGKLWVYYFDY"'''
    new_test = '''    fixture = ROOT / "data/regression/igfold_7nxb.pdb"
    receipt = json.loads(fixture.with_suffix(".receipt.json").read_text())
    assert digest(fixture) == receipt["sha256"]
    assert pdb_h3(fixture.read_text()) == "AKDGGKLWVYYFDY"'''
    assert old_test in text
    text = text.replace(old_test, new_test)
    target.write_text(text)
    target = DEST / "src/reconstruct_revision_references.py"
    text = target.read_text().replace("import io\n", "import io\nimport gzip\n")
    text = text.replace('(CACHE / leaf).read_bytes()', '(gzip.open(CACHE / "ab_split.csv.gz", "rb").read() if leaf == "ab_split.csv" else (CACHE / leaf).read_bytes())')
    text = text.replace('CACHE / "ab_split.csv", Path(__file__)', 'CACHE / "ab_split.csv.gz", Path(__file__)')
    begin = text.index('    report = f"""## Material Passport')
    end = text.index('## Result and scope', begin)
    text = text[:begin] + '    report = f"""# Reference reconstruction and alternate-atom policy\n\n' + text[end:]
    text = text.replace('metadata field, sequence, span, coordinate and null. No shared input, manuscript,\nresponse, checklist or saved model score was changed.',
                        'metadata field, sequence, span, coordinate and null. The frozen benchmark\nreferences and saved model scores are unchanged.')
    text = text.replace('UV_CACHE_DIR=/private/tmp/nile-reference-uv-cache uv run --offline --no-project --python 3.12 --with biopython==1.85 --with numpy==2.2.6',
                        'uv run --no-project --python 3.12 --with-requirements environments/reference.lock')
    begin = text.index('Local archive searches covered the project, Desktop, Downloads, /private/tmp and')
    end = text.index('\n"""', begin)
    text = text[:begin] + 'This reconstruction makes no claim to recover its historical predecessor and\nrecalculates no model scores. The upstream release and extracted member hashes\ndefine the reconstruction sources.' + text[end:]
    target.write_text(text)
    target = DEST / "src/recreate_reused_figures.py"
    text = target.read_text()
    begin, end = text.index('source = ROOT.parent'), text.index('assert source.is_file()')
    text = text[:begin] + 'source = ROOT / "paper_figure_sources/make_figures.py"\n' + text[end:]
    target.write_text(text)
    write_regression_fixture(source)

    for row in records:
        packaged = DEST / row["path"]
        row["packaged_sha256"] = digest(packaged)
        if row["transformation"] == "byte-identical" and row["packaged_sha256"] != row["source_sha256"]:
            row["transformation"] = "Public storage/portable-path adaptation; numerical algorithms retained; see docs/DEPENDENCIES.md"
    save(DEST / "provenance/source_selection.json", {"schema": 1, "files": records,
        "selection_basis": "Reviewed prepare_revision_handoff.py allowlist and dependency graph; no old verification clones copied",
        "excluded": ["private transcripts and ownership/response/checklist documents", "manuscripts and manuscript builders", "credentials", "model weights and native executables", "multi-GB archives, incomplete archive downloads and range parts", "old handoff trees and verification clones"]})
    paths = [p for p in (DEST / "outputs/revision_2026-10-16").rglob("*") if p.is_file()]
    paths += [DEST / "data/inputs/targets_vhh.json", DEST / "src/revision_training_audit.py"]
    paths += [DEST / p for p in ["outputs/structure/simplefold/metrics/per_target_scores.json", "outputs/structure/igfold/metrics/per_target_scores_igfold.json", "outputs/structure/nbb2/per_target_scores_nbb2.json"]]
    save(DEST / "outputs/revision_2026-10-17/overlap_extension/frozen_input_sha256.json", {str(p.relative_to(DEST)): digest(p) for p in sorted(paths)})
    save(DEST / "provenance/public_frozen_receipt_scope.json", {
        "receipt": "outputs/revision_2026-10-17/overlap_extension/frozen_input_sha256.json",
        "scope": "Every file in the curated public 16 October tree, frozen targets, public audit script and original three score JSONs",
        "historical_receipt_reused": False, "numeric_data_modified": False,
        "excluded_from_old_tree": ["private README.txt", "integrity_manifest.json referencing response documents", "raw ab_split.csv (losslessly gzipped)", "raw igfold_sabdab.zip (exact derived pool and checksum supplied)", "figure binaries (moved to figures/saved)"]})
    assert not (DEST / "inference/repeat_outputs").exists()
    print(json.dumps({"source_files": len(records), "embedding_shards": len(shards), "destination": str(DEST)}))


if __name__ == "__main__":
    main()
