"""TNP CDR3-conformation metrics (Therapeutic Nanobody Profiler, Commun Biol 9:344, 2026).

Definitions taken from the paper, not guessed:

  compactness = (number of residues in CDR3) / reach
  reach       = distance from the CDR3 anchor residues to the loop's furthest point
  anchors     = IMGT 105 and 117, i.e. the first and last CDR3 residues

  kink angles use backbone CA at IMGT 116, 117, 118:
      tau116   = bond angle     CA(115)-CA(116)-CA(117)
      alpha116 = dihedral       CA(115)-CA(116)-CA(117)-CA(118)
      kinked   : 0 < alpha116 < 120  and  85 < tau116 < 130
      extended : alpha116 < -100     and  100 < tau116 < 145

Index mapping: `cdr3 = [start, end)` spans IMGT 105..117 in hseq, so
IMGT 117 = end-1, 116 = end-2, 115 = end-3, 118 = end. IMGT insertions sit at
111/112, in the middle of the loop, so the C-terminal anchor arithmetic is safe.

TNP's own amber/red thresholds come from 36 clinical-stage nanobodies and live in a
Table-2 image we cannot read, so callers should percentile-calibrate on their own
population and say so.
"""
import numpy as np


def _ang(a, b, c):
    u, v = a - b, c - b
    cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(cos, -1, 1))))


def _dih(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def compactness(ca, span):
    """ca: (L,3) CA coords, may contain NaN. span: [start, end) into the chain.

    Returns (compactness, reach, n_res) or (None, None, n_res) if the loop is unusable.
    """
    lo, hi = span
    loop = ca[lo:hi]
    ok = np.isfinite(loop[:, 0])
    n_res = int(hi - lo)
    if ok.sum() < 3 or not (np.isfinite(ca[lo, 0]) and np.isfinite(ca[hi - 1, 0])):
        return None, None, n_res
    base = (ca[lo] + ca[hi - 1]) / 2.0           # anchors: IMGT 105 and 117
    reach = float(np.linalg.norm(loop[ok] - base, axis=-1).max())
    if reach <= 0:
        return None, None, n_res
    return n_res / reach, reach, n_res


def kink(ca, span):
    """Returns (alpha116, tau116, label) with label in {kinked, extended, other, na}."""
    hi = span[1]
    idx = [hi - 3, hi - 2, hi - 1, hi]            # IMGT 115, 116, 117, 118
    if idx[0] < 0 or idx[3] >= len(ca):
        return None, None, 'na'
    pts = ca[idx]
    if not np.isfinite(pts[:, 0]).all():
        return None, None, 'na'
    tau = _ang(pts[0], pts[1], pts[2])
    alpha = _dih(pts[0], pts[1], pts[2], pts[3])
    if 0 < alpha < 120 and 85 < tau < 130:
        lab = 'kinked'
    elif alpha < -100 and 100 < tau < 145:
        lab = 'extended'
    else:
        lab = 'other'
    return alpha, tau, lab


def band(value, lo_amber, hi_amber):
    """Percentile-calibrated flag band. TNP's absolute thresholds are unavailable."""
    if value is None:
        return 'na'
    if value < lo_amber:
        return 'amber_low'
    if value > hi_amber:
        return 'amber_high'
    return 'green'


def ca_from_hca(hca):
    """targets_vhh.json `hca` is [residue][N,CA,C,CB][xyz]; pull the CA column."""
    return np.array([[np.nan] * 3 if r is None or r[1] is None else r[1] for r in hca], float)


if __name__ == '__main__':  # self-check on synthetic geometry
    import math
    L, lo, hi = 20, 5, 15
    ca = np.zeros((L, 3))
    ca[:, 0] = np.arange(L) * 3.8                                  # straight chain
    c, r, n = compactness(ca, (lo, hi))
    assert n == 10 and abs(r - 4.5 * 3.8) < 1e-6, (r, n)           # 10 residues -> half-span 4.5 steps
    # a loop folded back on itself must have a smaller reach, hence larger compactness
    ca2 = ca.copy()
    ca2[lo:hi, 0] = ca[lo, 0] + np.abs(np.arange(hi - lo) - (hi - lo) / 2) * 1.0
    c2, r2, _ = compactness(ca2, (lo, hi))
    assert r2 < r and c2 > c, (r, r2, c, c2)
    assert band(0.5, 0.7, 1.3) == 'amber_low' and band(2.0, 0.7, 1.3) == 'amber_high'
    assert band(1.0, 0.7, 1.3) == 'green'
    print('tnp_metrics self-check OK')
