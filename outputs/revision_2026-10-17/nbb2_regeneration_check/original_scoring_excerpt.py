import biotite.structure.io.pdb as pdb_io

AA3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
          'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
          'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
REGIONS = ['cdr1','cdr2','cdr3','framework']

def kabsch(P, Q):
    Pm, Qm = P.mean(0), Q.mean(0)
    V, _, W = np.linalg.svd((P - Pm).T @ (Q - Qm))
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ W))])
    return V @ D @ W, Pm, Qm

rmsd = lambda a, b: float(np.sqrt(((a - b) ** 2).sum(-1).mean()))

rows = []
for t in TG:
    tid = t['pdb'] + t['H']
    try:
        f = f'{PRED}/{tid}.pdb'
        if not os.path.exists(f):
            R.setdefault('missing', []).append(tid); continue
        a = pdb_io.PDBFile.read(f).get_structure(model=1)
        a = a[np.isin(a.atom_name, ['N','CA','C','CB'])]

        seen, order = {}, []
        for ch, ri, ins, rn in zip(a.chain_id, a.res_id, a.ins_code, a.res_name):
            key = (str(ch), int(ri), str(ins))
            if key not in seen:
                seen[key] = len(order); order.append(AA3to1.get(str(rn), 'X'))
        pseq = ''.join(order)

        off = t['hseq'].find(pseq)
        if off < 0:
            R.setdefault('seq_mismatch', []).append(
                f"{tid}: pred {len(pseq)} vs ref {len(t['hseq'])}"); continue
        Ln = len(t['hseq'])
        pred = np.full((Ln, 4, 3), np.nan)
        ai = {'N':0,'CA':1,'C':2,'CB':3}
        for ch, ri, ins, an, xyz in zip(a.chain_id, a.res_id, a.ins_code, a.atom_name, a.coord):
            pred[seen[(str(ch), int(ri), str(ins))] + off, ai[str(an)]] = xyz
        ref = np.array([[[np.nan]*3 if v is None else v for v in r] for r in t['hca']], float)

        reg = {}
        for k in ('cdr1','cdr2','cdr3'):
            m = np.zeros(Ln, bool); lo, hi = t[k]; m[lo:hi] = True; reg[k] = m
        reg['framework'] = ~(reg['cdr1'] | reg['cdr2'] | reg['cdr3'])

        ok = np.isfinite(ref[..., 0]) & np.isfinite(pred[..., 0])
        if ok[:, 1].sum() < 60:
            R.setdefault('skipped', []).append(f'{tid}:coverage'); continue

        row = {'instance': t['instance'], 'tid': tid, 'L': Ln,
               'cdr3_len': int(t['cdr3'][1] - t['cdr3'][0]), 'n_atoms_matched': int(ok.sum()),
               'seq_offset': int(off), 'pred_len': len(pseq),
               'method': t['method'], 'holo': t['holo'], 'ab_split': t['ab_split'],
               'resolution': t['resolution']}
        for mode, fitmask in (('chain', np.ones(Ln, bool)), ('fw', reg['framework'])):
            sel = ok & fitmask[:, None]
            if sel.sum() < 60: continue
            Rot, Pm, Qm = kabsch(pred[sel], ref[sel])
            al = np.full_like(pred, np.nan)
            al[ok] = (pred[ok] - Pm) @ Rot + Qm
            for rn in REGIONS:
                s2 = ok & reg[rn][:, None]
                row[f'{mode}_{rn}'] = round(rmsd(al[s2], ref[s2]), 3) if s2.sum() >= 3 else None
        rows.append(row)
    except Exception as e:
        R.setdefault('score_err', []).append(f'{tid}:{type(e).__name__}:{str(e)[:90]}')

