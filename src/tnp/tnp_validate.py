"""Does the compactness bimodality correspond to the kinked/extended split?

Commun Biol 2023 (VHH CDR-H3 conformation vs VH germline): extended loops protrude into
solvent, kinked loops fold back onto framework. Folding back means a shorter reach, so
kinked should sit at HIGHER compactness. If that holds, the antimode of the compactness
histogram is a data-derived decision boundary, and we do not need TNP's unreadable Table 2.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
D = json.load(open(ROOT / 'data/inputs/tnp_crystals.json'))
rows = D['rows']
grp = {k: np.array([r['compact'] for r in rows if r['kink'] == k])
       for k in ('kinked', 'extended', 'other')}

print('compactness by kink class (crystals)')
for k, v in grp.items():
    print('  %-9s n=%3d  median %.3f  IQR %.3f-%.3f'
          % (k, len(v), np.median(v), np.percentile(v, 25), np.percentile(v, 75)))

# Mann-Whitney U via rank sum, normal approximation (no scipy on this machine)
a, b = grp['kinked'], grp['extended']
allv = np.concatenate([a, b])
r = np.argsort(np.argsort(allv)) + 1.0
U = r[:len(a)].sum() - len(a) * (len(a) + 1) / 2
mu = len(a) * len(b) / 2
sd = np.sqrt(len(a) * len(b) * (len(a) + len(b) + 1) / 12)
z = (U - mu) / sd
auc = U / (len(a) * len(b))
print('\nkinked vs extended:  AUC %.3f   z %.1f   (AUC 0.5 = no separation)' % (auc, z))

# antimode: the minimum of a smoothed histogram between the two peaks
h, e = np.histogram([r_['compact'] for r_ in rows], bins=40)
ctr = (e[:-1] + e[1:]) / 2
sm = np.convolve(h, np.ones(5) / 5, mode='same')
lo, hi = np.searchsorted(ctr, 0.85), np.searchsorted(ctr, 1.15)
anti = float(ctr[lo + int(np.argmin(sm[lo:hi]))])
print('antimode of compactness between the two peaks: %.3f' % anti)
side = lambda v: 'high' if v >= anti else 'low'
agree = sum(1 for x in rows if x['kink'] in ('kinked', 'extended')
            and ((x['kink'] == 'kinked') == (side(x['compact']) == 'high')))
n_cls = sum(1 for x in rows if x['kink'] in ('kinked', 'extended'))
print('antimode split agrees with the alpha/tau label on %d/%d = %.1f%%' % (agree, n_cls, 100 * agree / n_cls))

json.dump({'antimode': anti, 'p5': D['p5'], 'p95': D['p95'], 'auc_kink': auc},
          open(ROOT / 'data/inputs/tnp_boundary.json', 'w'), indent=1)
print('\nwrote tnp_boundary.json')
