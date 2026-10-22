from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import tomllib

E = Path('/home/dgamroth/workspaces/codex/modulation-storage-evidence')
A = Path('/home/dgamroth/workspaces/codex/AdaptiveOpticsSim-modulation-storage')
C = Path('/home/dgamroth/workspaces/codex/REVOLTCopperSim-modes')
D = A / 'benchmarks/results/platform/2026-10-21-modulation-storage'
R = D / 'raw'
R.mkdir(parents=True, exist_ok=True)
for label in ('baseline-final', 'changed'):
    source = E / label
    record = tomllib.loads((source/'validation.toml').read_text())
    assert record['complete']
    shutil.copytree(source, R/label, dirs_exist_ok=True)
    # Validate every archived scientific-array checksum, including phase tests.
    def verify(x):
        if isinstance(x, dict):
            if 'file' in x and 'sha256' in x:
                p = R/label/x['file']
                assert hashlib.sha256(p.read_bytes()).hexdigest() == x['sha256'], p
            for v in x.values(): verify(v)
        elif isinstance(x,list):
            for v in x: verify(v)
    verify(record)
for name in ('baseline-final.log','changed.log','cpu-tests.log','cuda-tests.log',
             'cpu-tests-final.log','copper-selector-tests.log','launcher-tests.log',
             'preparation-baseline-1440-256.toml','preparation-changed-1440-256.toml',
             'preparation-baseline-1440-256.log','preparation-changed-1440-256.log',
             'changed-initial-threshold.log','changed-initial-upload.log',
             'copper-modes.toml','copper-modes.log'):
    shutil.copy2(E/name,R/name)
for label in ('baseline','changed','modes'):
    dest=R/(label+'-env');dest.mkdir(exist_ok=True)
    for name in ('Project.toml','Manifest.toml'):
        shutil.copy2(E/(label+'-env')/name,dest/name)
shutil.copy2(E/'comparison.json',D/'comparison.json')
shutil.copy2(C/'benchmarks/copper_modes_cuda.jl',R/'copper_modes_cuda.jl')
shutil.copy2(Path('/home/dgamroth/workspaces/codex/REVOLTCopperSim.jl/benchmarks/copper_convergence.jl'),R/'copper_convergence.jl')
shutil.copy2(Path(__file__),R/'archive-results.py')
for name,command in (
    ('gpu.txt',['nvidia-smi','--query-gpu=name,driver_version,compute_cap,memory.total,power.limit','--format=csv']),
    ('system.txt',['uname','-a']),
    ('cpu.txt',['lscpu'])):
    (R/name).write_text(subprocess.check_output(command,text=True))
provenance={'baseline_aos_revision':'074f15ebeeb454a2e185aac58b2dc77a3a84894f',
            'implementation_commit':subprocess.check_output(['git','-C',str(A),'rev-parse','HEAD'],text=True).strip(),
            'launcher_commit':subprocess.check_output(['git','-C','/home/dgamroth/workspaces/codex/adaptive-optics-simulations-modes','rev-parse','HEAD'],text=True).strip(),
            'copper_modes_commit':subprocess.check_output(['git','-C',str(C),'rev-parse','HEAD'],text=True).strip(),
            'note':'Measurements preceded source commits; committed source content is the tested implementation. Source metadata and exact harness/environment checksums remain in raw TOML records.',
            'raw':{str(p.relative_to(D)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
                   for p in sorted(R.rglob('*')) if p.is_file()}}
(D/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
print('Archived',len(provenance['raw']),'raw artifacts')
