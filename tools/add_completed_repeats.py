"""Import completed coordinate-retaining repeats without changing primary data."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
RUNS = {'simplefold': 'simplefold_version1', 'nbb2': 'nbb2_version2'}


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()
    selection_path = ROOT / 'provenance/source_selection.json'
    selection = json.loads(selection_path.read_text())
    receipt = []
    for model, version in RUNS.items():
        original = source / 'outputs/revision_2026-10-23/repeat_runs' / model
        target = ROOT / 'inference/repeat_outputs' / version
        assert not target.exists(), 'Never overwrite a versioned repeat'
        assert json.loads((original / 'status.json').read_text())['status'] == 'COMPLETE'
        report = json.loads((original / 'verification.json').read_text())
        assert report['coordinate_scoring_exact_rows']
        assert not report['original_only'] and not report['repeat_only']
        provenance = json.loads((original / 'output/run_provenance.json').read_text())
        hashes = provenance['checkpoint_sha256']
        assert hashes and all(re.fullmatch('[0-9a-f]{64}', x) for x in hashes.values())
        if model == 'nbb2':
            assert {Path(p).name for p in hashes} == {f'nanobody_model_{i}' for i in range(1, 5)}
        paths = list((original / 'output').iterdir()) + list((original / 'source').iterdir())
        paths += [original / 'verification.json', original / 'status.json']
        for path in paths:
            assert path.is_file() and not path.is_symlink()
            relative = path.relative_to(original)
            destination = target / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            transformation = 'Byte-identical completed repeat artifact'
            if path.suffix == '.ipynb':
                nb = json.loads(path.read_text())
                for cell in nb['cells']:
                    if cell['cell_type'] == 'code':
                        cell['outputs'] = []
                        cell['execution_count'] = None
                destination.write_text(json.dumps(nb, indent=1) + '\n')
                transformation = 'Source cells preserved; notebook outputs and execution counts cleared'
            else:
                shutil.copyfile(path, destination)
            selection['files'].append({'path': str(destination.relative_to(ROOT)),
                'source': str(path.relative_to(source)), 'source_root': 'study',
                'source_sha256': sha(path), 'packaged_sha256': sha(destination),
                'transformation': transformation})
        name = 'per_target_scores.json' if model == 'simplefold' else 'per_target_scores_nbb2.json'
        original_scores = 'simplefold/metrics/' + name if model == 'simplefold' else 'nbb2/' + name
        primary = {r['tid']: r for r in json.loads((ROOT / 'outputs/structure' / original_scores).read_text())}
        repeated = {r['tid']: r for r in json.loads((target / 'output' / name).read_text())}
        assert set(primary) == set(repeated)
        with (target / 'comparison.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['tid', 'metric', 'primary_score', 'repeat_score', 'repeat_minus_primary'])
            for tid in sorted(primary):
                for metric in report['comparisons']:
                    a, b = primary[tid][metric], repeated[tid][metric]
                    writer.writerow([tid, metric, a, b, b-a])
        receipt.append({'run': version, 'coordinates': report['coordinate_files'],
            'archive_sha256': report['archive_sha256'], 'score_sha256': report['score_sha256'],
            'checkpoint_hashes_logged': len(hashes), 'primary_scores_replaced': False})
    selection_path.write_text(json.dumps(selection, indent=2) + '\n')
    (ROOT / 'provenance/completed_repeats.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
