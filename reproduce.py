"""Check benchmark files and reproduce saved analyses in an isolated working copy."""
from pathlib import Path, PurePosixPath
import argparse
import ast
import concurrent.futures
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import zipfile

ROOT = Path(__file__).resolve().parent


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def artifact_paths():
    ignored = {".git", ".venv", "__pycache__", ".ruff_cache", ".matplotlib"}
    return sorted(p for p in ROOT.rglob("*") if p.is_file() and not any(x in ignored for x in p.relative_to(ROOT).parts)
                  and p.name not in {"MANIFEST.sha256", "FILE_LIST.txt", ".DS_Store"}
                  and p != ROOT / "outputs/esm-probe/esm_probe_embeddings.npz")


def seal():
    paths = artifact_paths()
    (ROOT / "FILE_LIST.txt").write_text("\n".join(str(p.relative_to(ROOT)) for p in paths) + "\nFILE_LIST.txt\nMANIFEST.sha256\n")
    paths = sorted(paths + [ROOT / "FILE_LIST.txt"])
    (ROOT / "MANIFEST.sha256").write_text("".join(f"{digest(p)}  {p.relative_to(ROOT)}\n" for p in paths))
    print(json.dumps({"sealed_files": len(paths), "bytes": sum(p.stat().st_size for p in paths)}))


def verify_coordinate_archive(archive, entries):
    expected = {row["path"]: row for row in entries}
    assert entries and len(expected) == len(entries), "Empty or duplicate coordinate manifest"
    for name in expected:
        path = PurePosixPath(name)
        assert not path.is_absolute() and ".." not in path.parts and "\\" not in name, name
    with zipfile.ZipFile(archive) as stream:
        names = stream.namelist()
        assert len(names) == len(set(names)) and set(names) == set(expected), "ZIP member mismatch"
        assert stream.testzip() is None
        for name in names:
            payload = stream.read(name)
            assert len(payload) == expected[name]["bytes"], name
            assert hashlib.sha256(payload).hexdigest() == expected[name]["sha256"], name
    return len(entries)


def verify_completed_repeats():
    root = ROOT / "inference/repeat_outputs"
    if not root.exists():
        return []
    reports = []
    for run in sorted(p for p in root.iterdir() if p.is_dir()):
        assert json.loads((run / "status.json").read_text())["status"] == "COMPLETE", run.name
        files = run / "output"
        archives = list(files.glob("*_coordinates.zip"))
        assert len(archives) == 1, run.name
        archive = archives[0]
        model = archive.name.removesuffix("_coordinates.zip")
        assert model in {"nbb2", "simplefold"}
        entries = json.loads((files / "coordinate_manifest.json").read_text())
        count = verify_coordinate_archive(archive, entries)
        name = "per_target_scores_nbb2.json" if model == "nbb2" else "per_target_scores.json"
        saved = json.loads((files / name).read_text())
        tids = {row["tid"] for row in saved}
        assert len(tids) == len(saved) == count
        if model == "nbb2":
            assert tids == {PurePosixPath(row["path"]).stem for row in entries}
        report = json.loads((run / "verification.json").read_text())
        assert report["coordinate_files"] == count and report["archive_sha256"] == digest(archive)
        assert report["score_sha256"] == digest(files / name)
        assert report["coordinate_scoring_replay_max_absolute_difference"] <= 0.00100001
        provenance = json.loads((files / "run_provenance.json").read_text())
        assert provenance["input_sha256"] == digest(ROOT / "data/inputs/targets_vhh.json")
        hashes = provenance.get("checkpoint_sha256", {})
        assert isinstance(hashes, dict) and all(re.fullmatch(r"[0-9a-f]{64}", value) for value in hashes.values())
        primary_path = "nbb2/per_target_scores_nbb2.json" if model == "nbb2" else "simplefold/metrics/per_target_scores.json"
        primary = json.loads((ROOT / "outputs/structure" / primary_path).read_text())
        assert tids == {row["tid"] for row in primary} and count == (588 if model == "nbb2" else 597)
        assert report["coordinate_scoring_exact_rows"] and not report["original_only"] and not report["repeat_only"]
        if run.name != "nbb2_version1":
            assert hashes, "Completed newer run must include logged checkpoint hashes"
            if model == "nbb2":
                assert {PurePosixPath(p).name for p in hashes} == {f"nanobody_model_{i}" for i in range(1, 5)}
        reports.append({"run": run.name, "model": model, "coordinates": count,
                        "coordinate_member_hashes": "PASS", "checkpoint_hashes_logged": len(hashes),
                        "saved_coordinate_scoring_exact_rows": report["coordinate_scoring_exact_rows"]})
    assert reports, "Repeat directory contains no completed archives"
    return reports


def inventory():
    forbidden_parts = {"weights", "credentials", "verification_7p2dn81t", "omar_analysis_handoff", "word_qa", "supplied_project"}
    forbidden_names = {"integrate_revision_manuscript.py", "finalize_revision_publication.py", "prepare_revision_handoff.py", "kaggle.json"}
    files = artifact_paths()
    for path in files:
        rel = path.relative_to(ROOT)
        assert not path.is_symlink() and not forbidden_parts.intersection(rel.parts)
        assert path.name not in forbidden_names and path.suffix not in {".docx", ".part", ".bin", ".pt", ".pth"}
        assert path.stat().st_size < 100_000_000, f"Oversized public blob: {rel}"
        if path.suffix == ".py":
            ast.parse(path.read_text(), filename=str(rel))
        if path.suffix in {".py", ".json", ".jsonl", ".ipynb", ".md", ".txt", ".log", ".csv", ".yaml", ".toml"}:
            content = path.read_text()
            assert not re.search(r"hf_[A-Za-z0-9]{20,}|KGAT_[A-Za-z0-9_]{12,}|gh[pousr]_[A-Za-z0-9]{20,}|wandb_[A-Za-z0-9]{20,}", content), rel
            assert not re.search(r"/Users/[A-Za-z][A-Za-z0-9_-]*/", content), f"Private workstation path: {rel}"
        if path.suffix == ".ipynb":
            nb = json.loads(path.read_text())
            assert all(not c.get("outputs") and c.get("execution_count") is None for c in nb["cells"] if c["cell_type"] == "code")
    if (ROOT / "MANIFEST.sha256").exists():
        for line in (ROOT / "MANIFEST.sha256").read_text().splitlines():
            expected, name = line.split("  ", 1)
            assert digest(ROOT / name) == expected, name
        listed = {line.split("  ", 1)[1] for line in (ROOT / "MANIFEST.sha256").read_text().splitlines()}
        assert listed == {str(p.relative_to(ROOT)) for p in files} | {"FILE_LIST.txt"}
    frozen = json.loads((ROOT / "outputs/revision_2026-10-17/overlap_extension/frozen_input_sha256.json").read_text())
    assert all(digest(ROOT / name) == value for name, value in frozen.items())
    old = ROOT / "outputs/revision_2026-10-16"
    assert {str(p.relative_to(ROOT)) for p in old.rglob("*") if p.is_file()} <= set(frozen)
    targets = json.loads((ROOT / "data/inputs/targets_vhh.json").read_text())
    assert len(targets) == len({t["pdb"] + t["H"] for t in targets}) == 597
    for t in targets:
        for i in (1, 2, 3):
            assert t["hseq"][slice(*t[f"cdr{i}"])] == t[f"CDRH{i}"]
    reconstructed = json.loads((ROOT / "outputs/revision_2026-10-17/reference_reconstruction/targets_reconstructed_597.json").read_text())
    assert targets == reconstructed
    representatives = json.loads((ROOT / "outputs/revision_2026-10-17/reference_reconstruction/representatives_604.json").read_text())
    assert len(representatives) == 604
    for row in representatives:
        coordinate = ROOT / "outputs/revision_2026-10-17/reference_reconstruction/sources/coordinates" / (row["target"]["instance"] + ".cif")
        assert digest(coordinate) == row["evidence"]["source_coordinate_sha256"]
    score_files = ["simplefold/metrics/per_target_scores.json", "igfold/metrics/per_target_scores_igfold.json", "nbb2/per_target_scores_nbb2.json"]
    scores = [json.loads((ROOT / "outputs/structure" / p).read_text()) for p in score_files]
    assert [len(s) for s in scores] == [597, 597, 588]
    assert len(set.intersection(*({r["tid"] for r in s} for s in scores))) == 588
    repeats = verify_completed_repeats()
    return {"files": len(files), "targets": 597, "score_counts": [597, 597, 588], "common_targets": 588,
            "reference_json_exact": True, "public_frozen_entries": len(frozen), "configured_secret_scan": "PASS",
            "completed_repeat_archives": repeats}


def replay(args):
    reports = []
    out = ROOT / "verification"
    out.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="nanobody-public-verify-") as temp:
        work = Path(temp) / "repository"
        shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns(".git", ".venv", "__pycache__", ".matplotlib", "verification"))
        env = {**os.environ, "MPLCONFIGDIR": str(Path(temp) / "matplotlib"), "PYTHONUNBUFFERED": "1",
               "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}

        def run(label, script, environment="analysis", extra=(), comparisons=(), omit_keys=()):
            cmd = ["uv", "run", "--no-project", "--python", "3.13" if environment == "analysis" else "3.12",
                   "--with-requirements", f"environments/{environment}.lock", "python", script, *extra]
            if args.offline:
                cmd.insert(2, "--offline")
            started = time.monotonic()
            expected = {p: json.loads((work / p).read_text()) for p in comparisons}
            log = out / f"{label}.log"
            print(f"Running {label}", flush=True)
            with log.open("w") as stream:
                process = subprocess.run(cmd, cwd=work, stdout=stream, stderr=subprocess.STDOUT, env=env, timeout=3600)
            log_content = log.read_text().replace(str(work), "<fresh-repository>")
            if args.reference_archive:
                log_content = log_content.replace(str(args.reference_archive.resolve()), "<external-splits.tar.gz>")
            log.write_text(re.sub(r"/Users/[A-Za-z][A-Za-z0-9_-]*/", "<local-user>/", log_content))
            assert process.returncode == 0, f"{label} failed; see {log}"
            matches = {}
            geometry_comparison = None
            for relative, baseline in expected.items():
                actual = json.loads((work / relative).read_text())
                for key in omit_keys:
                    actual.pop(key, None)
                    baseline.pop(key, None)
                matches[relative] = actual == baseline
                if label == "crystal_labels":
                    assert {k: v for k, v in actual.items() if k != "rows"} == {k: v for k, v in baseline.items() if k != "rows"}
                    assert len(actual["rows"]) == len(baseline["rows"]) == 597
                    maximum = 0.0
                    differences = 0
                    for a, b in zip(actual["rows"], baseline["rows"]):
                        assert {k: v for k, v in a.items() if k not in {"alpha", "tau"}} == {k: v for k, v in b.items() if k not in {"alpha", "tau"}}
                        for key in ["alpha", "tau"]:
                            difference = abs(a[key] - b[key])
                            maximum = max(maximum, difference)
                            differences += int(difference != 0)
                    assert maximum <= 1e-12
                    geometry_comparison = {"labels_and_rounded_metrics_exact": True, "maximum_angle_difference_degrees": maximum,
                                           "angle_values_with_roundoff": differences, "absolute_angle_tolerance_degrees": 1e-12,
                                           "frozen_output_preserved": True}
                else:
                    assert matches[relative], f"{label} differs: {relative}"
            command_record = ["<external checksum-verified splits.tar.gz>" if args.reference_archive and c == str(args.reference_archive.resolve()) else c for c in cmd]
            result = {"name": label, "command": command_record, "exit_code": process.returncode,
                      "elapsed_seconds": round(time.monotonic() - started, 2), "exact_json_matches": matches,
                      "omitted_provenance_keys": list(omit_keys), "log": f"verification/{label}.log"}
            if geometry_comparison:
                result["geometry_comparison"] = geometry_comparison
            with (out / "completed_checks.jsonl").open("a") as stream:
                stream.write(json.dumps(result) + "\n")
            print(f"PASS {label}", flush=True)
            return result

        def probe():
            a = run("prepare_probe", "tools/prepare_probe_inputs.py")
            b = run("full_esm_probe", "src/reproduce_esm_probe.py", comparisons=("outputs/revision_2026-10-23/probe_reproduction/results.json",))
            report = json.loads((work / "outputs/revision_2026-10-23/probe_reproduction/verification.json").read_text())
            assert report["exact_saved_json_match"] and report["calls"] == 110
            assert report["packages"]["scikit-learn"] == "1.9.1"
            save = out / "probe_reproduction"
            save.mkdir(exist_ok=True)
            for filename in ["results.json", "verification.json", "out_of_fold_predictions.npz"]:
                shutil.copyfile(work / "outputs/revision_2026-10-23/probe_reproduction" / filename, save / filename)
            return [a, b]

        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            probe_job = pool.submit(probe) if not args.skip_probe and not args.reuse_probe else None
            if args.reuse_probe:
                previous = json.loads((out / "probe_reproduction/verification.json").read_text())
                baseline = json.loads((ROOT / "outputs/esm-probe/esm_probe_results.json").read_text())
                assert previous["exact_saved_json_match"] and previous["calls"] == 110
                assert json.loads((out / "probe_reproduction/results.json").read_text()) == baseline
                assert previous["input_sha256"]["src/esm_probe_analyze.py"] == digest(ROOT / "src/esm_probe_analyze.py")
                assert previous["input_sha256"]["data/inputs/targets_vhh.json"] == digest(ROOT / "data/inputs/targets_vhh.json")
                assert previous["input_sha256"]["outputs/esm-probe/esm_probe_meta.json"] == digest(ROOT / "outputs/esm-probe/esm_probe_meta.json")
                embedding_receipt = json.loads((ROOT / "data/esm_probe/embedding_manifest.json").read_text())
                assert previous["input_sha256"]["outputs/esm-probe/esm_probe_embeddings.npz"] == embedding_receipt["original_sha256"]
                assert previous["packages"]["scikit-learn"] == "1.9.1"
                reports.append({"name": "full_esm_probe", "status": "PASS", "execution": "Completed earlier in a fresh disposable copy during this assembly", "calls": 110,
                                "exact_saved_json_match": True, "receipt": "verification/probe_reproduction/verification.json"})
            reports.append(run("principal_analysis", "src/reviewer_revision_analysis.py", comparisons=("outputs/revision_2026-10-16/revision_statistics.json",)))
            reports.append(run("training_audit", "src/revision_training_audit.py", comparisons=("outputs/revision_2026-10-16/training_overlap_audit.json",)))
            reports.append(run("overlap_self_test", "src/revision_overlap_extension.py", "overlap", ("self-test",)))
            reports.append(run("overlap_summary", "src/revision_overlap_extension.py", "overlap", ("summarize",), ("outputs/revision_2026-10-17/overlap_extension/overlap_extension_summary.json",)))
            reports.append(run("simplefold_manifest", "src/revision_overlap_extension.py", "overlap", ("analyze-simplefold",), ("outputs/revision_2026-10-17/overlap_extension/simplefold_manifest_summary.json",)))
            reports.append(run("reference_self_check", "src/reconstruct_revision_references.py", "reference", ("--self-check",)))
            if args.reference_archive:
                reports.append(run("full_reference_reconstruction", "src/reconstruct_revision_references.py", "reference", ("--archive", str(args.reference_archive.resolve())), ("outputs/revision_2026-10-17/reference_reconstruction/targets_reconstructed_597.json", "outputs/revision_2026-10-17/reference_reconstruction/reference_policy_impact.json")))
                summary = json.loads((work / "outputs/revision_2026-10-17/reference_reconstruction/summary.json").read_text())
                assert summary["frozen_exact_matches"] == 597 and summary["archive_md5_verified"]
                (out / "reference_reconstruction.json").write_text(json.dumps(summary, indent=2) + "\n")
            reports.append(run("reference_sensitivity", "src/revision_reference_sensitivity.py", extra=("--expanded-cohort", "outputs/revision_2026-10-17/overlap_extension/extended_lower_risk_targets.csv"), comparisons=("outputs/revision_2026-10-17/reference_sensitivity/reference_sensitivity.json",), omit_keys=("reference_source_sha256",)))
            reports.append(run("reused_figures", "src/recreate_reused_figures.py", comparisons=("outputs/revision_2026-10-17/reused_figures/reproduction.json",)))
            reports.append(run("retained_structure_scoring", "src/verify_retained_structure_scoring.py"))
            reports.append(run("prepared_vhh_scoring", "src/score_prepared_vhh.py", extra=("--check",)))
            reports.append(run("complex_source_screen", "src/audit_complex_overlap.py", extra=("screen",), comparisons=("outputs/revision_2026-10-26/complex_overlap/antigen_matches.json",)))
            reports.append(run("complex_source_contacts", "src/audit_complex_overlap.py", extra=("contacts",), comparisons=("outputs/revision_2026-10-26/complex_overlap/summary.json", "outputs/revision_2026-10-26/complex_overlap/native_contacts.json", "outputs/revision_2026-10-26/complex_overlap/contact_extraction_errors.json")))
            reports.append(run("conformation_metrics", "src/tnp/tnp_metrics.py"))
            reports.append(run("crystal_labels", "src/tnp/tnp_crystals.py", comparisons=("data/inputs/tnp_crystals.json",)))
            reports.append(run("workflow_figure", "src/make_revision_workflow.py"))
            reports.append(run("archive_verifier_cli", "src/verify_archive_reruns.py", extra=("--help",)))
            if probe_job:
                reports.extend(probe_job.result())
        generated = list((work / "outputs/revision_2026-10-16").glob("benchmark_distributions_*"))
        generated += list((work / "outputs/revision_2026-10-17/reused_figures").glob("*.png"))
        generated += list((work / "outputs/revision_2026-10-17/reused_figures").glob("*.pdf"))
        generated += list((work / "outputs/revision_2026-10-23/publication").glob("workflow.*"))
        for path in generated:
            shutil.copyfile(path, out / path.name)
        result = {"status": "PASS", "fresh_disposable_copy": True, "source_study_modified": False,
                  "model_inference": False, "kaggle_kernels_launched": 0, "jobs": reports,
                  "full_reference_archive_replayed": bool(args.reference_archive),
                  "full_probe_repeated": not (args.skip_probe or args.reuse_probe),
                  "full_probe_reuse_validated": bool(args.reuse_probe)}
        (out / "reexecution.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seal", action="store_true")
    parser.add_argument("--replay", action="store_true")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--skip-probe", action="store_true")
    parser.add_argument("--reuse-probe", action="store_true", help="Reuse this assembly's previously completed and validated fresh-copy full probe")
    parser.add_argument("--reference-archive", type=Path)
    parser.add_argument("--optional-environment", action="store_true")
    args = parser.parse_args()
    assert not (args.skip_probe and args.reuse_probe)
    if args.optional_environment:
        import importlib.metadata
        import sys
        import biotite.structure.io.pdbx
        import scipy
        command = [sys.executable, "src/verify_archive_reruns.py", "--help"]
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=True)
        report = {"status": "PASS", "packages": {n: importlib.metadata.version(n) for n in ["numpy", "scipy", "biotite", "kaggle"]},
                  "biotite_pdbx_import": True, "archive_verifier_help": result.returncode == 0,
                  "remote_requests": 0, "kernels_launched": 0, "scope": "Optional stack import and CLI only; use the archive verifier for completed repeat data"}
        (ROOT / "verification/optional_environment.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report))
        return
    if args.seal:
        seal()
        return
    print(json.dumps(inventory()), flush=True)
    if args.replay:
        result = replay(args)
        print(json.dumps({"status": result["status"], "jobs": len(result["jobs"])}))


if __name__ == "__main__":
    main()
