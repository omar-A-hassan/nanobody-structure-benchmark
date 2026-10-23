"""Merge the SimpleFold and NanoBodyBuilder2 TNP arms.

The payoff question: SimpleFold is the only model here that reports a *distribution* of
compactness. Does that self-reported disagreement flag the cases where a TNP call is wrong --
including when the call was made by NanoBodyBuilder2, the model TNP actually uses?
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SF = {r['tid']: r for r in json.load(open(ROOT / 'outputs/tnp/sf/tnp_sf.json'))['rows']}
NB = {r['tid']: r for r in json.load(open(ROOT / 'outputs/tnp/specialists/tnp_specialists.json'))['rows']
      if r.get('nbb2_compact') is not None}
C = sorted(set(SF) & set(NB))
print('paired targets:', len(C))

rank = lambda v: np.argsort(np.argsort(np.asarray(v, float))).astype(float)
sp = lambda a, b: float(np.corrcoef(rank(a), rank(b))[0, 1])


def auc(pos, neg):
    a, b = np.asarray(pos, float), np.asarray(neg, float)
    if not len(a) or not len(b):
        return float('nan')
    r = rank(np.concatenate([a, b])) + 1
    return float((r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b)))


sf_flip = {t: SF[t]['sf_conf'] != SF[t]['cry_conf'] for t in C}
nb_flip = {t: NB[t]['nbb2_conf'] != NB[t]['cry_conf'] for t in C}
strad = {t: bool(SF[t]['straddles']) for t in C}

print('\nconformational call flips vs crystal')
print('  SimpleFold        %3d / %d = %.2f%%' % (sum(sf_flip.values()), len(C),
                                                 100 * sum(sf_flip.values()) / len(C)))
print('  NanoBodyBuilder2  %3d / %d = %.2f%%' % (sum(nb_flip.values()), len(C),
                                                 100 * sum(nb_flip.values()) / len(C)))
print('  both flip         %3d' % sum(1 for t in C if sf_flip[t] and nb_flip[t]))
print('  SimpleFold straddles (own 6 folds disagree)  %3d = %.2f%%'
      % (sum(strad.values()), 100 * sum(strad.values()) / len(C)))

print('\n--- 1. does the straddle flag predict SimpleFold\'s own wrong call? ---')
for lab, sel in (('straddles', True), ('agrees', False)):
    s = [t for t in C if strad[t] == sel]
    print('  %-10s n=%3d   flip rate %5.1f%%' % (lab, len(s), 100 * sum(sf_flip[t] for t in s) / len(s)))

print('\n--- 2. does it predict NanoBodyBuilder2\'s wrong call?  (cross-model) ---')
for lab, sel in (('straddles', True), ('agrees', False)):
    s = [t for t in C if strad[t] == sel]
    print('  %-10s n=%3d   NBB2 flip rate %5.1f%%' % (lab, len(s), 100 * sum(nb_flip[t] for t in s) / len(s)))

print('\n--- 3. continuous version: compactness spread as a score ---')
width = {t: SF[t]['compact_max'] - SF[t]['compact_min'] for t in C}
for nm, flip in (('SimpleFold', sf_flip), ('NanoBodyBuilder2', nb_flip)):
    a = auc([width[t] for t in C if flip[t]], [width[t] for t in C if not flip[t]])
    print('  AUC(compactness spread -> %-16s call is wrong) = %.3f' % (nm, a))

print('\n--- 4. control: is this just loop length? ---')
ln = {t: SF[t]['cdr3_len'] for t in C}
for nm, flip in (('SimpleFold', sf_flip), ('NanoBodyBuilder2', nb_flip)):
    print('  AUC(CDR3 length -> %-16s call is wrong)      = %.3f'
          % (nm, auc([ln[t] for t in C if flip[t]], [ln[t] for t in C if not flip[t]])))
print('  rho(compactness spread, CDR3 length) = %.3f' % sp([width[t] for t in C], [ln[t] for t in C]))

print('\n--- 5. does SimpleFold spread predict NBB2 compactness error? ---')
sprd = [SF[t]['spread_cdr3'] for t in C]
nberr = [abs(NB[t]['nbb2_dcompact']) for t in C]
sferr = [abs(SF[t]['sf_compact'] - SF[t]['cry_compact']) for t in C]
print('  rho(SF spread, SF  compactness error)   = %.3f' % sp(sprd, sferr))
print('  rho(SF spread, NBB2 compactness error)  = %.3f' % sp(sprd, nberr))
