"""Verify and unpack the TNP specialist result."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
R = json.load(open(ROOT / 'outputs/tnp/specialists/tnp_specialists.json'))
rows = [r for r in R['rows'] if r.get('nbb2_compact') is not None and r['cry_compact'] is not None]
print('n =', len(rows))

# 1. the identical-percentage red flag
fb = [r for r in rows if r['nbb2_band'] != r['cry_band']]
fc = [r for r in rows if r['nbb2_conf'] != r['cry_conf']]
print('\nband flips  %d (%.2f%%)' % (len(fb), 100 * len(fb) / len(rows)))
print('conf flips  %d (%.2f%%)' % (len(fc), 100 * len(fc) / len(rows)))
print('same targets? ', {r['tid'] for r in fb} == {r['tid'] for r in fc},
      '  overlap:', len({r['tid'] for r in fb} & {r['tid'] for r in fc}))

# 2. direction of the conformational flips
import collections
print('\nkinked/extended flip direction (crystal -> nbb2):')
for k, v in collections.Counter((r['cry_conf'], r['nbb2_conf']) for r in fc).most_common():
    print('   %-9s -> %-9s  %d' % (k[0], k[1], v))

# 3. how close to the boundary are the flippers?
anti = 0.961
d_all = np.array([abs(r['cry_compact'] - anti) for r in rows])
d_fl = np.array([abs(r['cry_compact'] - anti) for r in fc])
print('\ndistance of crystal compactness from the 0.961 antimode:')
print('   all targets   median %.4f' % np.median(d_all))
print('   flipped ones  median %.4f' % np.median(d_fl))

# 4. does structural error drive metric error?
rank = lambda v: np.argsort(np.argsort(np.asarray(v, float))).astype(float)
sp = lambda a, b: float(np.corrcoef(rank(a), rank(b))[0, 1])
have = [r for r in rows if r.get('nbb2_cdr3_rmsd') is not None]
rmsd = [r['nbb2_cdr3_rmsd'] for r in have]
dcomp = [abs(r['nbb2_dcompact']) for r in have]
print('\nrho(CDR3 RMSD, |dcompactness|) = %.3f   n=%d' % (sp(rmsd, dcomp), len(have)))
flip_rmsd = [r['nbb2_cdr3_rmsd'] for r in have if r['nbb2_conf'] != r['cry_conf']]
keep_rmsd = [r['nbb2_cdr3_rmsd'] for r in have if r['nbb2_conf'] == r['cry_conf']]
print('   median CDR3 RMSD, flipped  %.2f A  (n=%d)' % (np.median(flip_rmsd), len(flip_rmsd)))
print('   median CDR3 RMSD, kept     %.2f A  (n=%d)' % (np.median(keep_rmsd), len(keep_rmsd)))

# 5. flip rate as a function of structural accuracy
print('\nflip rate by CDR3 RMSD bin:')
for lo, hi in ((0, 1.5), (1.5, 2.5), (2.5, 4), (4, 6), (6, 99)):
    s = [r for r in have if lo <= r['nbb2_cdr3_rmsd'] < hi]
    if not s:
        continue
    f = sum(1 for r in s if r['nbb2_conf'] != r['cry_conf'])
    print('   %4.1f-%-4.1f A  n=%3d  flips %3d  = %5.1f%%' % (lo, hi, len(s), f, 100 * f / len(s)))
