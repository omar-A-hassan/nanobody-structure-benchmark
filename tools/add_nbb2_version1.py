"""Append the completed NBB2 v1 snapshot; never regenerate notebooks or launch jobs."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "inference/repeat_outputs/nbb2_version1"
ORIGIN = "outputs/revision_2026-10-23/repeat_runs/nbb2_version1"


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--finalize", action="store_true")
    args = parser.parse_args()
    source = args.source.resolve()
    selection_path = ROOT / "provenance/source_selection.json"
    selection = json.loads(selection_path.read_text())
    if not args.finalize:
        assert not RUN.exists(), "Never overwrite an existing versioned repeat"
        snapshot = ROOT / "provenance/core_snapshot"
        assert not snapshot.exists()
        snapshot.mkdir()
        for name in ["MANIFEST.sha256", "FILE_LIST.txt"]:
            shutil.copyfile(ROOT / name, snapshot / name)
        protected = {}
        for line in (snapshot / "MANIFEST.sha256").read_text().splitlines():
            digest, name = line.split("  ", 1)
            if name.startswith(("verification/", "data/", "outputs/", "figures/", "paper_figure_sources/", "environments/")) or name in {"pyproject.toml", "uv.lock"}:
                assert sha(ROOT / name) == digest, name
                protected[name] = digest
        save(snapshot / "protected_sha256.json", protected)
        assert json.loads((source / ORIGIN / "status.json").read_text())["status"] == "COMPLETE"
        required = ["repeat.ipynb", "nbb2-benchmark-archive-20261023.ipynb", "kernel-metadata.json", "status.json", "verification.json"]
        required += ["output/" + name for name in ["coordinate_manifest.json", "environment.txt", "nbb2-benchmark-archive-20261023.log",
                    "nbb2_coordinates.zip", "per_target_scores_nbb2.csv", "per_target_scores_nbb2.json", "results_nbb2.json",
                    "run_provenance.json", "worker0.log", "worker1.log"]]
        copies = [(f"{ORIGIN}/{name}", f"inference/repeat_outputs/nbb2_version1/{name}") for name in required]
        copies += [("src/verify_archive_reruns.py", "src/verify_archive_reruns.py"),
                   ("src/verify_archive_reruns.py", "provenance/original_scripts/verify_archive_reruns.py"),
                   ("src/prepare_archive_reruns.py", "src/prepare_archive_reruns.py")]
        for original, relative in copies:
            path, target = source / original, ROOT / relative
            assert path.is_file() and not path.is_symlink()
            if path.suffix in {".py", ".ipynb", ".json", ".csv", ".txt", ".log"}:
                assert not re.search(r"hf_[A-Za-z0-9]{20,}|KGAT_[A-Za-z0-9_]{12,}|gh[pousr]_[A-Za-z0-9]{20,}|wandb_[A-Za-z0-9]{20,}", path.read_text()), original
                assert not re.search(r"/Users/[A-Za-z][A-Za-z0-9_-]*/", path.read_text()), original
            target.parent.mkdir(parents=True, exist_ok=True)
            note = "byte-identical completed v1 artifact or current source"
            if path.suffix == ".ipynb":
                nb = json.loads(path.read_text())
                for cell in nb["cells"]:
                    if cell["cell_type"] == "code":
                        cell["outputs"] = []
                        cell["execution_count"] = None
                target.write_text(json.dumps(nb, indent=1) + "\n")
                note = "Notebook source cells preserved; outputs/execution counts cleared"
            else:
                shutil.copyfile(path, target)
            row = {"path": relative, "source": original, "source_root": "study", "source_sha256": sha(path),
                   "packaged_sha256": sha(target), "transformation": note}
            existing = next((r for r in selection["files"] if r["path"] == relative), None)
            if existing is None:
                selection["files"].append(row)
            else:
                existing.update(row)
        save(selection_path, selection)
        print(json.dumps({"copied_files": len(copies), "uncompressed_coordinate_duplicates": 0}))
        return

    protected = json.loads((ROOT / "provenance/core_snapshot/protected_sha256.json").read_text())
    assert all(sha(ROOT / name) == value for name, value in protected.items()), "Core artifacts changed"
    report = json.loads((ROOT / "verification/nbb2_version1_reexecution.json").read_text())
    original = json.loads((RUN / "verification.json").read_text())
    for key in ["coordinate_files", "archive_sha256", "score_sha256", "common_targets", "original_only", "repeat_only",
                "coordinate_scoring_replay_max_absolute_difference", "coordinate_scoring_exact_rows", "score_errors", "comparisons"]:
        assert report[key] == original[key], key
    assert report["coordinate_files"] == report["common_targets"] == 588
    assert report["coordinate_scoring_exact_rows"] and report["coordinate_scoring_replay_max_absolute_difference"] == 0
    provenance = json.loads((RUN / "output/run_provenance.json").read_text())
    assert provenance["checkpoint_sha256"] == {}
    assert report["checkpoint_provenance"]["status"] == "not_recorded"
    legacy = ROOT / "src/consolidate_stats.py"
    if legacy.exists():
        destination = ROOT / "provenance/original_scripts/consolidate_stats.py"
        assert not destination.exists()
        legacy.rename(destination)
        row = next(r for r in selection["files"] if r["path"] == "src/consolidate_stats.py")
        row["path"] = str(destination.relative_to(ROOT))
        row["transformation"] = "Byte-identical legacy aggregator relocated to provenance; not an active analysis entry point"
    save(RUN / "comparison.json", {"run": "nbb2_version1", "primary_replaced": False, "common_targets": 588,
         "original_only": report["original_only"], "repeat_only": report["repeat_only"], "comparisons": report["comparisons"],
         "coordinate_scoring_exact_rows": True, "checkpoint_hashes_known": False})
    primary = {r["tid"]: r for r in json.loads((ROOT / "outputs/structure/nbb2/per_target_scores_nbb2.json").read_text())}
    repeat = {r["tid"]: r for r in json.loads((RUN / "output/per_target_scores_nbb2.json").read_text())}
    assert set(primary) == set(repeat) and len(primary) == 588
    with (RUN / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["tid", "metric", "primary_score", "repeat_score", "repeat_minus_primary"])
        writer.writeheader()
        for tid in sorted(primary):
            for metric in report["comparisons"]:
                a, b = primary[tid][metric], repeat[tid][metric]
                writer.writerow({"tid": tid, "metric": metric, "primary_score": a, "repeat_score": b, "repeat_minus_primary": b - a})
    for row in selection["files"]:
        if row["path"] == "src/verify_archive_reruns.py":
            row["packaged_sha256"] = sha(ROOT / row["path"])
            row["transformation"] = "Public versioned local-run/report options, disposable coordinate extraction, and explicit checkpoint logging validation; original scorer unchanged"
    save(selection_path, selection)
    copied = [r for r in selection["files"] if r["source"].startswith(ORIGIN + "/")]
    for row in copied:
        assert sha(ROOT / row["path"]) == row["packaged_sha256"]
        assert sha(source / row["source"]) == row["source_sha256"]
    command = [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"]
    tests = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    (ROOT / "verification/nbb2_version1_boundary_tests.log").write_text(tests.stdout + tests.stderr)
    assert tests.returncode == 0, "Repeat boundary tests failed"
    save(ROOT / "provenance/nbb2_version1_extension.json", {"status": "PASS", "run": "nbb2_version1",
         "source": ORIGIN, "files_copied": len(copied), "coordinate_archive_sha256": report["archive_sha256"],
         "coordinate_count": 588, "same_primary_ids": True, "coordinate_scoring_exact_rows": True,
         "checkpoint_hashes_known": False, "checkpoint_logging_gap": "V1 extension filter omitted extensionless nanobody_model_* weights; raw checkpoint_sha256 remains empty",
         "primary_replaced": False, "core_protected_files_unchanged": len(protected), "prior_core_check_count": 15,
         "prior_core_check_report_sha256": sha(ROOT / "verification/reexecution.json"),
         "local_reexecution_report": "verification/nbb2_version1_reexecution.json", "kernels_launched": 0,
         "new_model_inference": False, "v2_and_simplefold_results_added": False,
         "boundary_tests": {"exit_code": tests.returncode, "log": "verification/nbb2_version1_boundary_tests.log"},
         "legacy_aggregator": "provenance/original_scripts/consolidate_stats.py"})
    print(json.dumps({"status": "PASS", "core_files_preserved": len(protected), "repeat_coordinates": 588}))


if __name__ == "__main__":
    main()
