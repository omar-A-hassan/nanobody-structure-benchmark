"""TNP CDR3 compactness on the 597 SAbDab2 crystal references."""
import json
from pathlib import Path
import numpy as np
from tnp_metrics import compactness, kink, ca_from_hca

ROOT = Path(__file__).resolve().parents[2]
T = json.load(open(ROOT / 'data/inputs/targets_vhh.json'))
rows, lab = [], {'kinked': 0, 'extended': 0, 'other': 0, 'na': 0}
for t in T:
    ca = ca_from_hca(t['hca'])
    c, r, n = compactness(ca, t['cdr3'])
    if c is None:
        continue
    a, tau, l = kink(ca, t['cdr3'])
    lab[l] += 1
    rows.append({'tid': f"{t['pdb']}{t['H']}", 'compact': round(c, 4),
                 'reach': round(r, 3), 'n_res': n, 'alpha': a, 'tau': tau, 'kink': l})

C = np.array([r['compact'] for r in rows])
R = np.array([r['reach'] for r in rows])
N = np.array([r['n_res'] for r in rows])
rank = lambda v: np.argsort(np.argsort(v))
sp = lambda a, b: float(np.corrcoef(rank(a), rank(b))[0, 1])

print(f'597 crystals -> {len(rows)} scored')
print('compactness  median %.3f   p5 %.3f   p95 %.3f   range %.3f-%.3f'
      % (np.median(C), np.percentile(C, 5), np.percentile(C, 95), C.min(), C.max()))
print('reach (A)    median %.2f   range %.2f-%.2f' % (np.median(R), R.min(), R.max()))
print('CDR3 length  median %d   range %d-%d' % (np.median(N), N.min(), N.max()))
print('rho(compact, length) %.3f    rho(compact, reach) %.3f' % (sp(C, N), sp(C, R)))
print('kink rule:', lab, '-> classified %.1f%%' % (100 * (lab['kinked'] + lab['extended']) / len(rows)))

h, e = np.histogram(C, bins=22)
print('\ncompactness histogram (is it bimodal?)')
for k in range(22):
    print('  %.2f-%.2f | %s %d' % (e[k], e[k + 1], '#' * max(0, h[k] // 2), h[k]))

json.dump({'p5': float(np.percentile(C, 5)), 'p95': float(np.percentile(C, 95)),
           'median': float(np.median(C)), 'rows': rows},
          open(ROOT / 'data/inputs/tnp_crystals.json', 'w'), indent=1)
print('\nwrote tnp_crystals.json')
