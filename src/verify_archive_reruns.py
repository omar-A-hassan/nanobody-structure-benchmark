"""Verify retained repeat coordinates and compare repeat scores with the benchmark."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import glob
import hashlib
import importlib.metadata as metadata
import json
import os
import re
import tempfile
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/revision_2026-10-23/repeat_runs'
REFS = {
    'simplefold': 'potatoguy44/simplefold-benchmark-archive-20261023',
    'nbb2': 'potatoguy44/nbb2-benchmark-archive-20261023',
}


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def code_cells(path):
    return [''.join(c['source']) for c in json.loads(path.read_text())['cells']
            if c['cell_type'] == 'code']


def checkpoint_provenance(hashes, model, required=False):
    assert isinstance(hashes, dict)
    assert all(isinstance(value, str) and re.fullmatch(r'[0-9a-f]{64}', value) for value in hashes.values())
    numbered = {Path(name).name for name in hashes if Path(name).name.startswith('nanobody_model_')}
    if required:
        assert hashes, 'Checkpoint hashes were not logged for this run'
        if model == 'nbb2':
            assert numbered == {f'nanobody_model_{i}' for i in range(1, 5)}, 'Four NBB2 weight hashes are required'
    return {'status': 'logged' if hashes else 'not_recorded', 'hashes_logged': len(hashes),
            'nbb_weight_hashes_logged': len(numbered), 'required': required,
            'weights_redistributed': False, 'weights_independently_rehashed': False}


def verify(model, folder=None, report_path=None, require_checkpoint_hashes=False, require_exact=False):
    folder = folder or OUT / model
    with tempfile.TemporaryDirectory(prefix='nanobody-repeat-coordinates-') as temporary:
        return verify_coordinates(model, folder, Path(temporary), report_path, require_checkpoint_hashes, require_exact)


def verify_coordinates(model, folder, target_root, report_path, require_checkpoint_hashes, require_exact):
    files = folder / 'output'
    if (folder / 'status.json').exists():
        assert json.loads((folder / 'status.json').read_text())['status'] == 'COMPLETE'
    provenance = json.loads((files / 'run_provenance.json').read_text())
    assert provenance['input_sha256'] == sha(ROOT / 'data/inputs/targets_vhh.json')
    checkpoint = checkpoint_provenance(provenance.get('checkpoint_sha256', {}), model, require_checkpoint_hashes)
    archive = files / (model + '_coordinates.zip')
    entries = json.loads((files / 'coordinate_manifest.json').read_text())
    expected = {e['path']: e for e in entries}
    assert entries and len(expected) == len(entries), 'Duplicate/empty coordinate manifest'
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert len(z.namelist()) == len(set(z.namelist())), 'Duplicate ZIP member'
        assert set(z.namelist()) == set(expected)
        for name in z.namelist():
            relative = Path(name)
            assert not relative.is_absolute() and '..' not in relative.parts and '\\' not in name
            data = z.read(name)
            entry = expected[name]
            assert len(data) == entry['bytes']
            assert hashlib.sha256(data).hexdigest() == entry['sha256']
            destination = target_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
    targets = json.loads((ROOT / 'data/inputs/targets_vhh.json').read_text())
    original_name = 'sf-nbench-full' if model == 'simplefold' else 'nbb2-bench'
    source = next((ROOT / 'outputs/revision_2026-10-23/recovery' / original_name / 'source').glob('*.ipynb'))
    scoring = next(c for c in code_cells(source) if 'def kabsch' in c)
    # Execute the original coordinate-to-score implementation, without its writes.
    prefix = scoring.split("R['rows'] = rows")[0]
    namespace = {'np': np, 'os': os, 'glob': glob, 'TG': targets, 'R': {},
                 'SCRATCH': str(target_root), 'PRED': str(target_root),
                 'tid_of': {t['instance']: t['pdb'] + t['H'] for t in targets}}
    exec(compile(prefix, str(source) + ':scoring', 'exec'), namespace)
    recomputed = namespace['rows']
    score_name = 'per_target_scores.json' if model == 'simplefold' else 'per_target_scores_nbb2.json'
    saved = json.loads((files / score_name).read_text())
    original = json.loads((ROOT / 'outputs/revision_2026-10-23/recovery' / original_name / 'output' / score_name).read_text())
    fields = [m + '_' + r for m in ['chain', 'fw']
              for r in ['cdr1', 'cdr2', 'cdr3', 'framework']]
    saved_by = {r['tid']: r for r in saved}
    old_by = {r['tid']: r for r in original}
    replay_by = {r['tid']: r for r in recomputed}
    assert len(saved_by) == len(saved) == len(replay_by) == len(recomputed) == len(entries)
    assert set(saved_by) == set(replay_by)
    replay_differences = [abs(replay_by[t][f] - saved_by[t][f])
                          for t in saved_by for f in fields]
    assert max(replay_differences, default=0) <= 0.00100001
    if require_exact:
        assert saved == recomputed, 'Full saved coordinate-to-score rows differ'
    shared = sorted(set(saved_by) & set(old_by))
    comparisons = {}
    for f in fields:
        new = np.array([saved_by[t][f] for t in shared])
        old = np.array([old_by[t][f] for t in shared])
        diff = new - old
        comparisons[f] = {'original_median': float(np.median(old)),
                          'repeat_median': float(np.median(new)),
                          'median_paired_change': float(np.median(diff)),
                          'median_absolute_change': float(np.median(abs(diff))),
                          'maximum_absolute_change': float(np.max(abs(diff))),
                          'exact_score_matches': int(np.sum(diff == 0))}
    report = {'checked_utc': datetime.now(timezone.utc).isoformat(), 'model': model,
              'run_type': 'Documented repeat; primary scores remain unchanged',
              'coordinate_files': len(entries), 'archive_sha256': sha(archive),
              'score_sha256': sha(files / score_name), 'common_targets': len(shared),
              'original_only': sorted(set(old_by) - set(saved_by)),
              'repeat_only': sorted(set(saved_by) - set(old_by)),
              'coordinate_scoring_replay_max_absolute_difference': max(replay_differences, default=0),
              'coordinate_scoring_exact_rows': saved == recomputed,
              'score_errors': namespace['R'], 'comparisons': comparisons,
              'historical_coordinate_identity_established': False,
              'checkpoint_provenance': checkpoint,
              'primary_scores_replaced': False,
              'original_scoring_notebook_sha256': sha(source),
              'packages': {name: metadata.version(name) for name in ['numpy', 'scipy', 'biotite']}}
    report_path = report_path or folder / 'verification.json'
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--model', choices=list(REFS))
    parser.add_argument('--run-directory', type=Path, help='Completed local versioned run; never downloaded into')
    parser.add_argument('--report', type=Path, help='Separate local verification JSON; required with --run-directory')
    parser.add_argument('--require-checkpoint-hashes', action='store_true', help='Require complete hash logging for a newer run; omit for NBB2 v1')
    parser.add_argument('--require-exact', action='store_true', help='Require exact saved/replayed rows, beyond score-rounding tolerance')
    args = parser.parse_args()
    if args.run_directory:
        if not args.model or not args.report or args.download:
            parser.error('--run-directory requires --model and --report, and excludes --download')
        assert args.report.resolve() != (args.run_directory / 'verification.json').resolve(), 'Preserve the imported run verification'
        verify(args.model, args.run_directory.resolve(), args.report, args.require_checkpoint_hashes, args.require_exact)
        return
    api = None
    if args.download:
        from kaggle.api.kaggle_api_extended import KaggleApi
        api = KaggleApi()
        api.authenticate()
    for model in [args.model] if args.model else REFS:
        folder = OUT / model
        folder.mkdir(parents=True, exist_ok=True)
        if api:
            status = json.loads(str(api.kernels_status(REFS[model])))
            (folder / 'status.json').write_text(json.dumps(status, indent=2) + '\n')
            print(model, status['status'], flush=True)
            if status['status'] != 'COMPLETE':
                continue
            api.kernels_pull(REFS[model],path=str(folder / 'source'),metadata=True,quiet=True)
            api.kernels_output(REFS[model], str(folder / 'output'),
                               file_pattern=r'\.(json|csv|txt|log|zip|tar)$', quiet=True)
        verify(model, require_checkpoint_hashes=args.require_checkpoint_hashes, require_exact=args.require_exact)


if __name__ == '__main__':
    main()
