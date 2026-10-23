"""Build repeat-run notebooks that retain coordinates and execution records."""
from pathlib import Path
import ast
import argparse
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "kaggle/kernels/revision-archive-20261023"
OUT.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT / "outputs/revision_2026-10-23/recovery"


def replace_once(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)


provenance = '''
import hashlib, importlib.metadata, platform
def sha_file(path):
    with open(path, 'rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()
R['run_type'] = 'documented repeat; original saved scores are not replaced'
R['python'] = sys.version
R['platform'] = platform.platform()
R['input_sha256'] = sha_file(DATASET if 'DATASET' in globals() else _ds[0])
R['packages'] = {d.metadata['Name']: d.version for d in importlib.metadata.distributions()}
open(f'{WORK}/environment.txt', 'w').write(sh(f'{sys.executable} -m pip freeze').stdout)
json.dump(R, open(f'{WORK}/run_provenance.json', 'w'), indent=2, default=str)
'''

archive_sf = '''
import zipfile, hashlib
coordinate_files = sorted(glob.glob(f'{SCRATCH}/out*/**/*.cif', recursive=True))
assert len(coordinate_files) == len(TG), (len(coordinate_files), len(TG))
coordinate_manifest = []
with zipfile.ZipFile(f'{WORK}/simplefold_coordinates.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
    for p in coordinate_files:
        name = os.path.relpath(p, SCRATCH)
        archive.write(p, name)
        coordinate_manifest.append({'path': name, 'sha256': sha_file(p), 'bytes': os.path.getsize(p)})
json.dump(coordinate_manifest, open(f'{WORK}/coordinate_manifest.json', 'w'), indent=2)
R['checkpoint_sha256'] = {os.path.relpath(p, SCRATCH): sha_file(p) for p in glob.glob(f'{SCRATCH}/artifacts/**/*', recursive=True) if os.path.isfile(p)}
json.dump(R, open(f'{WORK}/run_provenance.json', 'w'), indent=2, default=str)
note('coordinates retained')
'''

archive_nb = '''
import zipfile, hashlib
coordinate_files = sorted(glob.glob(f'{PRED}/*.pdb'))
assert len(coordinate_files) >= 580, len(coordinate_files)
coordinate_manifest = []
with zipfile.ZipFile(f'{WORK}/nbb2_coordinates.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
    for p in coordinate_files:
        name = os.path.basename(p)
        archive.write(p, name)
        coordinate_manifest.append({'path': name, 'sha256': sha_file(p), 'bytes': os.path.getsize(p)})
json.dump(coordinate_manifest, open(f'{WORK}/coordinate_manifest.json', 'w'), indent=2)
from ImmuneBuilder import NanoBodyBuilder2
import inspect
package_dir = os.path.dirname(inspect.getfile(NanoBodyBuilder2))
R['source_sha256'] = {os.path.relpath(p, package_dir): sha_file(p) for p in glob.glob(package_dir+'/**/*.py', recursive=True)}
R['checkpoint_sha256'] = {os.path.relpath(p, package_dir): sha_file(p) for p in glob.glob(package_dir+'/**/*', recursive=True) if os.path.isfile(p) and (p.endswith(('.pt','.pth')) or os.path.basename(p).startswith('nanobody_model_'))}
assert len([p for p in R['checkpoint_sha256'] if os.path.basename(p).startswith('nanobody_model_')]) == 4
json.dump(R, open(f'{WORK}/run_provenance.json', 'w'), indent=2, default=str)
note('coordinates retained')
'''

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--model',choices=['simplefold','nbb2'])
args=parser.parse_args()
sha = None
if args.model != 'nbb2':
    sha = subprocess.check_output(['git', 'ls-remote', 'https://github.com/omar-A-hassan/ml-simplefold', 'refs/heads/fix-multichain-fasta-parsing'], text=True).split()[0]
    assert len(sha) == 40
for model, original in [('simplefold', 'sf-nbench-full'), ('nbb2', 'nbb2-bench')]:
    if args.model and model != args.model: continue
    source = next((SOURCE / original / 'source').glob('*.ipynb'))
    notebook = json.loads(source.read_text())
    codes = [c for c in notebook['cells'] if c['cell_type'] == 'code']
    for c in notebook['cells']:
        if c['cell_type'] == 'code':
            c['source'] = ''.join(c['source'])
            c['outputs'] = []
            c['execution_count'] = None
        else:
            c['source'] = ['Repeat monomer benchmark with retained coordinates and execution records.\n']
    if model == 'simplefold':
        codes[1]['source'] = replace_once(codes[1]['source'],
            "must(f'{sys.executable} -m pip install -q -e {SCRATCH}/repo', 'pip install simplefold')",
            f"must(f'git -C {{SCRATCH}}/repo checkout -q {sha}', 'pin source')\n" +
            "must(f'{sys.executable} -m pip install -q -e {SCRATCH}/repo', 'pip install simplefold')")
        codes[1]['source'] += provenance + '''
R['source_commit'] = must(f'git -C {SCRATCH}/repo rev-parse HEAD', 'source revision').stdout.strip()
must(f'git -C {SCRATCH}/repo archive -o {WORK}/simplefold_source.tar HEAD', 'archive source')
'''
        codes[5]['source'] = replace_once(codes[5]['source'], "t0 = time.time()\nprocs = []", "R['worker_fasta_order'] = {str(g): glob.glob(f'{SCRATCH}/sh{g}/*.fasta') for g in (0, 1)}\nt0 = time.time()\nprocs = []")
        codes[5]['source'] = replace_once(codes[5]['source'], "out = p.communicate()[0]", "out = p.communicate()[0]\n    open(f'{WORK}/worker{g}.log', 'w').write(out)")
        codes[8]['source'] = archive_sf
    else:
        codes[1]['source'] = replace_once(codes[1]['source'], 'install -q ImmuneBuilder', 'install -q ImmuneBuilder==1.2')
        codes[1]['source'] += provenance
        codes[2]['source'] = replace_once(codes[2]['source'], "out = p.communicate()[0]", "out = p.communicate()[0]\n    open(f'{WORK}/worker{g}.log', 'w').write(out)")
        codes[3]['source'] = replace_once(codes[3]['source'], 'shutil.rmtree(PRED, ignore_errors=True)', archive_nb)
    for i,c in enumerate(notebook['cells']):
        c['id'] = f'repeat-{model}-{i:03d}'
        if c['cell_type'] == 'code':
            ast.parse(c['source'])
    folder = OUT / model
    folder.mkdir(exist_ok=True)
    (folder / 'repeat.ipynb').write_text(json.dumps(notebook, indent=1)+'\n')
    metadata = {'id': f'potatoguy44/{model}-benchmark-archive-20261023', 'title': f'{model} benchmark archive 20261023',
        'code_file': 'repeat.ipynb', 'language': 'python', 'kernel_type': 'notebook', 'is_private': True,
        'enable_gpu': True, 'machine_shape': 'NvidiaTeslaT4', 'enable_internet': True,
        'dataset_sources': ['potatoguy44/sabdab2-nanobodies']}
    (folder/'kernel-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (folder/'source_receipt.json').write_text(json.dumps({'source_notebook_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'simplefold_repeat_source_commit': sha if model=='simplefold' else None,
        'scope': 'Coordinate retention and execution provenance added; original regional scoring retained'},indent=2)+'\n')
    print(folder)
