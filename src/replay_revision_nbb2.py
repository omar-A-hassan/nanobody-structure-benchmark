"""Three-target CPU inference-to-score check; never replaces historical scores."""
from pathlib import Path
import hashlib
import importlib.metadata as metadata
import inspect
import json
import os
import platform
import shutil
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/revision_2026-10-17/nbb2_regeneration_check"
OUT.mkdir(parents=True, exist_ok=True)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


target_path = ROOT / "data/inputs/targets_vhh.json"
score_path = ROOT / "outputs/structure/nbb2/per_target_scores_nbb2.json"
notebook_path = ROOT / "kaggle/kernels/nbb2-bench/nbb2-bench.ipynb"
original = {r["tid"]: r for r in json.loads(score_path.read_text())}
targets = json.loads(target_path.read_text())
selected = sorted((t for t in targets if str(t["holo"]).lower() == "false"
                   and t["pdb"] + t["H"] in original), key=lambda t: t["pdb"] + t["H"])[:3]
assert [t["pdb"] + t["H"] for t in selected] == ["7BNWB", "7DJXA", "7DJYB"]
assert all(t["hseq"][t["cdr3"][0]:t["cdr3"][1]] == t["CDRH3"] for t in selected)
frozen = {str(p.relative_to(ROOT)): digest(p) for p in (target_path, score_path, notebook_path)}
save("inputs.json", selected)
report = {"scope": "New, isolated three-target standalone VHH CPU regeneration; not antigen-context inference or docking",
          "selection": "Lexicographically first three source-free VHHs with retained NBB2 scores; no error-based selection",
          "historical_prediction_recovery": False, "original_inputs_sha256": frozen,
          "python": sys.version, "platform": platform.platform(), "argv": sys.argv,
          "script_sha256": digest(Path(__file__)),
          "status": "STARTED", "targets": []}
started = time.monotonic()
try:
    import numpy as np
    import torch
    import hmmer
    from ImmuneBuilder import NanoBodyBuilder2
    torch.set_num_threads(2)
    torch.manual_seed(12345)
    packages = {d.metadata["Name"]: d.version for d in metadata.distributions()}
    report["packages"] = packages
    (OUT / "requirements-nbb-replay.lock").write_text("\n".join(f"{name}=={version}" for name, version in sorted(packages.items(), key=lambda x: x[0].lower())) + "\n")
    # The wheel's launcher renames its binary, breaking SearchIO's program-name
    # detection. Use the unchanged packaged native binary with its real name.
    packaged = Path(hmmer.__file__).parent / "data/bin/hmmscan"
    native = OUT / "bin/hmmscan"
    native.parent.mkdir(exist_ok=True)
    shutil.copyfile(packaged, native)
    native.chmod(0o755)
    assert digest(packaged) == digest(native)
    os.environ["PATH"] = str(native.parent) + os.pathsep + os.environ.get("PATH", "")
    assert shutil.which("hmmscan") == str(native)
    probe = subprocess.run(["hmmscan", "-h"], text=True, capture_output=True, check=True)
    assert "HMMER 3.4" in probe.stdout and probe.stdout.startswith("# hmmscan ::")
    report["hmmer_banner"] = probe.stdout.splitlines()[:4]
    report["hmmer_native_sha256"] = digest(native)
    from ImmuneBuilder.sequence_checks import number_sequences
    report["numbering_preflight"] = []
    for target in selected:
        numbered = number_sequences({"H": target["hseq"]}, allowed_species=None, scheme="imgt")["H"]
        h3 = "".join(aa for (number, _), aa in numbered if 105 <= number <= 117)
        assert h3 == target["CDRH3"]
        report["numbering_preflight"].append({"tid": target["pdb"] + target["H"], "h3_matches_frozen": True})
    weights = OUT / "weights"
    predictions = OUT / "predictions"
    weights.mkdir(exist_ok=True)
    predictions.mkdir(exist_ok=True)
    nb = NanoBodyBuilder2(weights_dir=str(weights), numbering_scheme="imgt")
    assert nb.device == "cpu", "This check is deliberately CPU-only"
    package_dir = Path(inspect.getfile(NanoBodyBuilder2)).parent
    report["model_source_sha256"] = {str(p.relative_to(package_dir)): digest(p) for p in sorted(package_dir.rglob("*.py"))}
    report["weight_sha256"] = {p.name: digest(p) for p in sorted(weights.iterdir()) if p.is_file()}
    for target in selected:
        tid = target["pdb"] + target["H"]
        path = predictions / (tid + ".pdb")
        assert not path.exists(), f"Do not silently reuse or overwrite a previous run: {path}"
        t0 = time.monotonic()
        nb.predict({"H": target["hseq"]}).save(str(path))
        assert path.is_file() and path.stat().st_size > 0
        report["targets"].append({"tid": tid, "seconds": round(time.monotonic() - t0, 2),
                                  "prediction_sha256": digest(path), "prediction": str(path.relative_to(ROOT))})
        save("run_receipt.json", report)
        print("predicted", tid, report["targets"][-1]["seconds"], "seconds", flush=True)

    notebook = json.loads(notebook_path.read_text())
    code = "\n".join(c["source"] if isinstance(c["source"], str) else "".join(c["source"])
                     for c in notebook["cells"] if c["cell_type"] == "code")
    start = "import biotite.structure.io.pdb as pdb_io"
    end = "R['rows'] = rows"
    assert code.count(start) == code.count(end) == 1
    excerpt = code.split(start, 1)[1].split(end, 1)[0]
    excerpt = start + excerpt
    assert "rmtree" not in excerpt and "json.dump" not in excerpt
    (OUT / "original_scoring_excerpt.py").write_text(excerpt)
    env = {"np": np, "os": os, "TG": selected, "PRED": str(predictions), "R": {}}
    exec(compile(excerpt, str(notebook_path) + ":original_scoring_excerpt", "exec"), env)
    rows = env["rows"]
    assert len(rows) == 3 and not env["R"], env["R"]
    keys = [f"{mode}_{region}" for mode in ("chain", "fw") for region in ("cdr1", "cdr2", "cdr3", "framework")]
    comparisons = [{"tid": r["tid"], "original": {k: original[r["tid"]][k] for k in keys},
                    "new": {k: r[k] for k in keys},
                    "new_minus_original": {k: round(r[k] - original[r["tid"]][k], 3) for k in keys}}
                   for r in rows]
    save("new_scores.json", rows)
    save("comparison.json", comparisons)
    report["maximum_absolute_regional_score_difference"] = max(abs(r["new_minus_original"][k]) for r in comparisons for k in keys)
    report["scorer_excerpt_sha256"] = digest(OUT / "original_scoring_excerpt.py")
    report["status"] = "PASS_NEW_PIPELINE_CHECK_NOT_HISTORICAL_REPRODUCTION"
except Exception:
    report["status"] = "FAILED"
    report["error"] = traceback.format_exc()
    raise
finally:
    report["elapsed_seconds"] = round(time.monotonic() - started, 2)
    assert all(digest(ROOT / p) == h for p, h in frozen.items()), "Historical inputs changed"
    save("run_receipt.json", report)
    print(json.dumps({k: report[k] for k in ("status", "elapsed_seconds", "targets")}), flush=True)
