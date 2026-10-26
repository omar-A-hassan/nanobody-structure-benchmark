"""Repeat the stored ESM2 probe and retain out-of-fold predictions."""
from pathlib import Path
import ast
import hashlib
import importlib.metadata
import json
import os

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
OUT = ROOT / 'outputs/revision_2026-10-23/probe_reproduction'
OUT.mkdir(parents=True, exist_ok=True)
source_path = ROOT / 'src/esm_probe_analyze.py'
source = source_path.read_text()
anchor = '# ---- exploratory thresholds'
assert source.count(anchor) == 1
prefix = source.split(anchor)[0]
calls = []


def record_prediction(original, estimator, X, y, **kwargs):
    prediction = original(estimator, X, y, **kwargs)
    calls.append({'prediction': prediction, 'labels': np.asarray(y),
                  'model': repr(estimator), 'groups': np.asarray(kwargs['groups'])})
    return prediction


injection = '''
_original_cross_val_predict = cross_val_predict
def cross_val_predict(estimator, X, y, **kwargs):
    return record_prediction(_original_cross_val_predict, estimator, X, y, **kwargs)
'''
prefix = prefix.replace('RNG = np.random.default_rng(0)', injection+'\nRNG = np.random.default_rng(0)')
scope = {'record_prediction': record_prediction}
exec(compile(prefix, str(source_path), 'exec'), scope)
scope['linear_pass'] = bool(scope['best']['auc'] >= .70 and scope['gap'] >= .15)
scope['mlp_pass'] = bool(scope['auc_mlp'] >= .70 and scope['gap_mlp'] >= .15)
tree = ast.parse(source)
assignment = next(n for n in tree.body if isinstance(n,ast.Assign)
                  and any(isinstance(t,ast.Name) and t.id=='out' for t in n.targets))
exec(compile(ast.Module(body=[assignment],type_ignores=[]),str(source_path),'exec'),scope)
actual = scope['out']
baseline = json.loads((ROOT/'outputs/esm-probe/esm_probe_results.json').read_text())
np.savez_compressed(OUT/'out_of_fold_predictions.npz', **{
    f'{i}_{field}': row[field] for i,row in enumerate(calls)
    for field in ['prediction','labels','groups']})
(OUT/'results.json').write_text(json.dumps(actual,indent=2)+'\n')
inputs = ['src/esm_probe_analyze.py','outputs/esm-probe/esm_probe_embeddings.npz',
          'outputs/esm-probe/esm_probe_meta.json','data/inputs/targets_vhh.json']
report = {'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
          'packages':{n:importlib.metadata.version(n) for n in ['numpy','scipy','scikit-learn']},
          'calls':len(calls),'exact_saved_json_match':actual==baseline,
          'different_fields':{k:{'saved':baseline.get(k),'repeat':v} for k,v in actual.items() if v!=baseline.get(k)},
          'best_layer':actual['best_layer'],'best_layer_auc':actual['best_layer_auc'],
          'scope':'Repeated fits from the retained embeddings; no model inference or new target population',
          'models':[c['model'] for c in calls]}
(OUT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='models'},indent=2))
