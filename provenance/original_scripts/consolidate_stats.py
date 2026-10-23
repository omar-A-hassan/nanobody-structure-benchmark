"""Re-derive every number cited in docs/ from the output files. Nothing from memory."""
import json
import random
import statistics as st

import numpy as np

OUT = {}
L = lambda p: json.load(open(p))
rows_of = lambda d: d if isinstance(d, list) else d['rows']


def rank(v):
    return np.argsort(np.argsort(np.asarray(v, float))).astype(float)


def sp(a, b):
    return float(np.corrcoef(rank(a), rank(b))[0, 1])


def boot_med(d, n=4000, seed=1):
    random.seed(seed)
    v = sorted(st.median(random.choices(d, k=len(d))) for _ in range(n))
    return round(st.median(d), 3), (round(v[int(.025 * n)], 3), round(v[int(.975 * n)], 3))


# ---------------------------------------------------------------- 1. benchmark
sf = {r['instance']: r for r in rows_of(L('outputs/structure/simplefold/metrics/per_target_scores.json'))}
nb = {r['instance']: r for r in rows_of(L('outputs/structure/nbb2/per_target_scores_nbb2.json'))}
ig = {r['instance']: r for r in rows_of(L('outputs/structure/igfold/metrics/per_target_scores_igfold.json'))}
C = sorted(set(sf) & set(nb) & set(ig))
K = 'chain_cdr3'
b = OUT['benchmark'] = {'n_sf': len(sf), 'n_nbb2': len(nb), 'n_igfold': len(ig), 'n_paired': len(C)}

b['medians_all597_sf'] = {k: round(st.median([r[k] for r in sf.values()]), 3)
                          for k in ('chain_cdr1', 'chain_cdr2', 'chain_cdr3', 'chain_framework',
                                    'fw_cdr1', 'fw_cdr2', 'fw_cdr3', 'fw_framework')}
b['paired_medians'] = {n: {k: round(st.median([d[i][k] for i in C]), 3)
                           for k in ('chain_cdr3', 'chain_framework')}
                       for n, d in (('simplefold', sf), ('igfold', ig), ('nbb2', nb))}
for n, d in (('igfold', ig), ('nbb2', nb)):
    diff = [sf[i][K] - d[i][K] for i in C]
    m, ci = boot_med(diff)
    b[f'paired_diff_sf_minus_{n}'] = {'median': m, 'ci95': list(ci),
                                      'sf_wins_pct': round(100 * sum(1 for x in diff if x < 0) / len(diff), 1)}
b['pct_over_6A'] = {n: round(100 * sum(1 for i in C if d[i][K] > 6) / len(C), 1)
                    for n, d in (('simplefold', sf), ('igfold', ig), ('nbb2', nb))}
b['percentiles_cdr3'] = {n: {f'p{q}': round(float(np.percentile([d[i][K] for i in C], q)), 2)
                             for q in (25, 50, 75, 90, 95)} | {'max': round(max(d[i][K] for i in C), 2)}
                         for n, d in (('simplefold', sf), ('igfold', ig), ('nbb2', nb))}
fails = {n: {i for i in C if d[i][K] > 6} for n, d in (('simplefold', sf), ('igfold', ig), ('nbb2', nb))}
b['failure_overlap'] = {'sf_only': len(fails['simplefold'] - fails['igfold'] - fails['nbb2']),
                        'all_three': len(fails['simplefold'] & fails['igfold'] & fails['nbb2']),
                        'both_specialists': len(fails['igfold'] & fails['nbb2'])}
# Pearson is inflated here by two targets whose deposited chain is stretched (9G13A, 8V9XA),
# so report all three views rather than the flattering one.
_a = np.array([sf[i]['chain_framework'] for i in C])
_b = np.array([nb[i]['chain_framework'] for i in C])
_m = (_a < 5) & (_b < 5)
b['framework_err_corr_sf_nbb2'] = {
    'pearson_all': round(float(np.corrcoef(_a, _b)[0, 1]), 3),
    'pearson_excl_2_extreme': round(float(np.corrcoef(_a[_m], _b[_m])[0, 1]), 3),
    'spearman_all': round(sp(_a, _b), 3),
    'n_excluded': int((~_m).sum())}
# model-independent stratification (by NBB2 error)
b['strata_by_nbb2'] = {}
for lab, pred in (('nbb2_solves_lt2.5', lambda i: nb[i][K] < 2.5),
                  ('nbb2_middling_2.5_6', lambda i: 2.5 <= nb[i][K] < 6),
                  ('nbb2_struggles_gt6', lambda i: nb[i][K] >= 6)):
    S = [i for i in C if pred(i)]
    b['strata_by_nbb2'][lab] = {
        'n': len(S),
        'simplefold': round(st.median([sf[i][K] for i in S]), 3),
        'igfold': round(st.median([ig[i][K] for i in S]), 3),
        'nbb2': round(st.median([nb[i][K] for i in S]), 3),
        'sf_beats_nbb2_pct': round(100 * sum(1 for i in S if sf[i][K] < nb[i][K]) / len(S), 1),
        'sf_beats_igfold_pct': round(100 * sum(1 for i in S if sf[i][K] < ig[i][K]) / len(S), 1)}

# ---------------------------------------------------------------- 2. hit rate
h = L('outputs/hitrate/run/hitrate_results.json')
OUT['hitrate'] = {'elapsed_min': h['elapsed_min'], 'n_samples': h['n_samples'], 'taus': h['taus'],
                  'model': h['model'], 'seed': h['seed'],
                  'arms': {t: {'hard': {k: a['hard'][k] for k in
                                        ('n', 'hit@1.5', 'hit@2.0', 'hit@2.5', 'hit@3.0',
                                         'fw_best_med', 'fw_med_med', 'bad_pep', 'clash')},
                               'easy': {k: a['easy'][k] for k in ('n', 'hit@2.0', 'fw_best_med', 'bad_pep')},
                               'best_of_k': a['hard']['best_of_k']}
                           for t, a in h['arms'].items()},
                  'verdict': h['verdict']}

# ---------------------------------------------------------------- 3. TNP
cry = L('data/inputs/tnp_crystals.json')
OUT['tnp_boundary'] = L('data/inputs/tnp_boundary.json')
OUT['tnp_crystals'] = {'p5': round(cry['p5'], 3), 'p95': round(cry['p95'], 3),
                       'median': round(cry['median'], 3), 'n': len(cry['rows'])}
spec = L('outputs/tnp/specialists/tnp_specialists.json')
sfj = L('outputs/tnp/sf/tnp_sf.json')
OUT['tnp'] = {'nbb2_summary': spec['summary']['nbb2'],
              'nbb2_built': spec['nbb2_built'], 'nbb2_refused': spec['nbb2_refused'],
              'igfold_built': spec['igfold_built'],
              'sf': sfj['tnp'], 'sf_elapsed_min': sfj['elapsed_min'],
              'sf_chunk_rc': sfj['chunk_rc']}
OUT['replication'] = sfj['replication']

# ---------------------------------------------------------------- 4. trajectory / TFM
OUT['tfm'] = L('outputs/tabular/sf-tfm-only-597/tfm_results.json')
tr = L('outputs/trajectory/full-597/traj_results.json')
trr = tr['rows']
OUT['trajectory_597'] = {'n': len(trr), 'tau': tr['tau'], 'steps': tr['steps'],
                         'n_samples': tr['n_samples'], 'log_every': tr['log_every'],
                         'rho_vs_err': {f: round(sp([r[f] for r in trr], [r['err'] for r in trr]), 3)
                                        for f in ('spread_cdr3', 'adiv_final', 'adiv_half',
                                                  'adiv_q1', 'div_final', 'cdr3_len')},
                         'median_err': round(st.median([r['err'] for r in trr]), 3),
                         'median_err_oracle': round(st.median([r['err_oracle'] for r in trr]), 3),
                         'median_err_mean': round(st.median([r['err_mean'] for r in trr]), 3)}

# ---------------------------------------------------------------- 5. ensemble / hard
e = L('outputs/ensemble/run/ensemble_results.json')
OUT['ensemble_run'] = {k: e[k] for k in e if k not in ('rows',)}
hd = L('outputs/hard/run/hard_results.json')
hrows = hd['rows']
OUT['hard_run'] = {'n': len(hrows),
                   'spread_median': {g: round(st.median([r['spread_cdr3'] for r in hrows if r['group'] == g]), 3)
                                     for g in ('hard', 'easy')},
                   'err_median': {g: round(st.median([r['err'] for r in hrows if r['group'] == g]), 2)
                                  for g in ('hard', 'easy')}}

json.dump(OUT, open('docs/consolidated_stats.json', 'w'), indent=1, default=str)
print(json.dumps(OUT, indent=1, default=str))
