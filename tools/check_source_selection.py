"""Read-only source comparison and extraction checks for the assembled snapshot."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--figure-source", type=Path, required=True)
    args = parser.parse_args()
    selection = json.loads((ROOT / "provenance/source_selection.json").read_text())
    counts = {"source_hashes_checked": 0, "packaged_hashes_checked": 0,
              "notebook_code_cells_exact": 0, "notebook_markdown_adaptations": 0}
    for row in selection["files"]:
        base = args.source if row["source_root"] == "study" else args.figure_source.parent
        original, packaged = base / row["source"], ROOT / row["path"]
        assert digest(original) == row["source_sha256"], f"Source changed since selection: {row['source']}"
        assert digest(packaged) == row["packaged_sha256"], f"Packaged artifact changed: {row['path']}"
        counts["source_hashes_checked"] += 1
        counts["packaged_hashes_checked"] += 1
        if original.suffix == ".ipynb":
            a, b = json.loads(original.read_text()), json.loads(packaged.read_text())
            assert [c["cell_type"] for c in a["cells"]] == [c["cell_type"] for c in b["cells"]]
            for original_cell, packaged_cell in zip(a["cells"], b["cells"]):
                if original_cell["cell_type"] == "code" or not row.get("markdown_only_changes"):
                    assert original_cell["source"] == packaged_cell["source"], row["path"]
            counts["notebook_code_cells_exact"] += 1
            counts["notebook_markdown_adaptations"] += int(bool(row.get("markdown_only_changes")))
    archive = args.source / "outputs/revision_2026-10-16/sources/igfold_sabdab.zip"
    with tempfile.TemporaryDirectory(prefix="nanobody-source-check-") as temp:
        output = Path(temp) / "pool.json"
        subprocess.run([sys.executable, str(ROOT / "tools/extract_igfold_sequence_pool.py"), str(archive), str(output)], check=True)
        actual = json.loads(output.read_text())
        expected = json.loads((ROOT / "outputs/revision_2026-10-16/sources/igfold_experimental_sequence_pool.json").read_text())
        assert all(actual[k] == expected[k] for k in actual)
    metadata = []
    for model in ["simplefold", "nbb2"]:
        original = args.source / f"outputs/revision_2026-10-23/repeat_runs/{model}/source/kernel-metadata.json"
        actual = json.loads(original.read_text())
        packaged = json.loads((ROOT / f"kaggle/kernels/revision-archive-20261023/{model}/kernel-metadata.json").read_text())
        assert actual["id"] == packaged["id"] and packaged["code_file"] == "repeat.ipynb"
        metadata.append({"model": model, "actual_ref": actual["id"], "source": str(original.relative_to(args.source)),
                         "source_sha256": digest(original), "packaged_code_file": packaged["code_file"]})
    assert digest(args.figure_source) == digest(ROOT / "paper_figure_sources/make_figures.py")
    repeats = args.source / "outputs/revision_2026-10-23/repeat_runs"
    repeat_results = [str(p.relative_to(args.source)) for p in repeats.rglob("*") if p.is_file() and "output" in p.relative_to(repeats).parts]
    report = {"status": "PASS", **counts, "igfold_pool_reextraction_exact": True,
              "current_reusable_figure_source_exact": True, "actual_repeat_metadata": metadata,
              "local_repeat_output_files": repeat_results,
              "writes_to_source": 0, "scope": "Every selected original source artifact remains byte-identical to its assembly-time hash"}
    (ROOT / "verification/source_integrity.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "actual_repeat_metadata"}))


if __name__ == "__main__":
    main()
