"""Summarize retained repeat predictions without replacing primary scores."""
from pathlib import Path
import hashlib
import json
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
FIELDS = [f'{fit}_{region}' for fit in ['chain', 'fw']
          for region in ['cdr1', 'cdr2', 'cdr3', 'framework']]
PRIMARY = {'simplefold': 'simplefold/metrics/per_target_scores.json',
           'igfold': 'igfold/metrics/per_target_scores_igfold.json',
           'nbb2': 'nbb2/per_target_scores_nbb2.json'}


def load(path):
    return json.loads(path.read_text())


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def compare(primary, repeated, ids):
    result = {}
    for field in FIELDS:
        a = [primary[t][field] for t in ids]
        b = [repeated[t][field] for t in ids]
        changes = [y-x for x, y in zip(a, b)]
        result[field] = {'primary_median': median(a), 'repeat_median': median(b),
            'median_paired_change': median(changes),
            'median_absolute_change': median(abs(d) for d in changes),
            'maximum_absolute_change': max(abs(d) for d in changes),
            'exact_values': sum(d == 0 for d in changes), 'n': len(ids)}
    return result


def main():
    primary = {m: {r['tid']: r for r in load(ROOT / 'outputs/structure' / p)}
               for m, p in PRIMARY.items()}
    common = sorted(set.intersection(*(set(rows) for rows in primary.values())))
    assert len(common) == 588
    runs = {}
    for model, version in [('simplefold', 1), ('nbb2', 1), ('nbb2', 2)]:
        name = f'{model}_version{version}'
        folder = ROOT / 'inference/repeat_outputs' / name
        report = load(ROOT / 'verification' / (name + '_reexecution.json'))
        assert report['coordinate_scoring_exact_rows']
        assert report['coordinate_scoring_replay_max_absolute_difference'] == 0
        path = folder / 'output' / Path(PRIMARY[model]).name
        repeated = {r['tid']: r for r in load(path)}
        assert set(repeated) == set(primary[model])
        ids = sorted(repeated)
        provenance = load(folder / 'output/run_provenance.json')
        runs[name] = {'model': model, 'version': version, 'retained_coordinates': len(ids),
            'all_model_targets': compare(primary[model], repeated, ids),
            'common_588': compare(primary[model], repeated, common),
            'coordinate_score_replay_exact': True,
            'checkpoint_hashes_logged': len(provenance['checkpoint_sha256']),
            'primary_score_sha256': sha(ROOT / 'outputs/structure' / PRIMARY[model]),
            'repeat_score_sha256': sha(path), 'archive_sha256': report['archive_sha256'],
            'stage_seconds': provenance.get('fold_sec', provenance.get('build_sec')),
            'primary_scores_replaced': False}
    output = ROOT / 'outputs/revision_2026-10-23/repeat_summary.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({'common_targets': common, 'runs': runs}, indent=2) + '\n')
    print(json.dumps({k: {'n': v['retained_coordinates'],
        'whole_chain_h3': v['common_588']['chain_cdr3']} for k, v in runs.items()}, indent=2))


if __name__ == '__main__':
    main()
