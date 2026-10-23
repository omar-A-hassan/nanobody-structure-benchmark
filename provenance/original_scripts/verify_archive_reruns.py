"""Verify retained repeat coordinates and compare repeat scores with the benchmark."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import glob
import hashlib
import json
import os
import re
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


def verify(model):
    folder = OUT / model
    files = folder / 'output'
    archive = files / (model + '_coordinates.zip')
    entries = json.loads((files / 'coordinate_manifest.json').read_text())
    target_root = folder / 'coordinates'
    target_root.mkdir(exist_ok=True)
    expected = {e['path']: e for e in entries}
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert set(z.namelist()) == set(expected)
        for name in z.namelist():
            relative = Path(name)
            assert not relative.is_absolute() and '..' not in relative.parts
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
    assert set(saved_by) == set(replay_by)
    replay_differences = [abs(replay_by[t][f] - saved_by[t][f])
                          for t in saved_by for f in fields]
    assert max(replay_differences, default=0) <= 0.00100001
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
              'historical_coordinate_identity_established': False}
    (folder / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--model', choices=list(REFS))
    args = parser.parse_args()
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
        verify(model)


if __name__ == '__main__':
    main()
