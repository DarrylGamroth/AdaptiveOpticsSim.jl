"""Summarize raw completed-node observations and fixed-input storage diagnostics.
Usage: python summarize_pupil_reuse.py EVIDENCE OUTPUT.json
"""
import hashlib
import json
from pathlib import Path
import sys
import tomllib
import numpy as np


def record(path):
    return tomllib.loads(path.read_text())


def rate(root, meta):
    path=root/meta['rate_file']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==meta['output_sha256']
    return np.fromfile(path,dtype='<f4').astype(np.float64)


def summarize(root):
    result={'timing_scope':'Unpaced completed optical node; shared host/uncontrolled clocks; small deltas are descriptive',
            'performance':{},'storage_diagnostics':{}}
    for case in ('copper_fast','gpi_default','copper_default','gpi_fine','copper_fine'):
        before,after=record(root/f'before-{case}.toml'),record(root/f'after-{case}.toml')
        data={'before_ms':before['median_ms'],'after_ms':after['median_ms'],
              'dense_phase_bytes':before['phase_bytes'],'factor_bytes':after['phase_bytes'],
              'pupil_phase_cache_bytes':after['pupil_phase_bytes'],
              'observed_pool_before_bytes':before['device_used_bytes'],
              'observed_pool_after_bytes':after['device_used_bytes']}
        for phase in ('before','after'):
            path=root/f'repeat-{phase}-{case}.toml'
            if path.exists():data[f'repeat_{phase}_ms']=record(path)['median_ms']
        result['performance'][case]=data
    for case in ('copper_fast','gpi_fine','copper_fine'):
        dense=record(root/f'dense-{case}.toml');reference=rate(root,dense)
        data={}
        for name,prefix in [('dense','dense'),('separable','after'),('computed','computed')]:
            meta=record(root/f'{prefix}-{case}.toml');values=rate(root,meta)
            relative=float(np.linalg.norm(values-reference)/np.linalg.norm(reference))
            assert relative<2e-6 and all(x==0 for x in meta['host_allocated_bytes'])
            data[name]={'median_ms':meta['median_ms'],'storage_bytes':meta['phase_bytes'],
                        'relative_l2':relative,'identical_to_dense':bool(np.array_equal(values,reference))}
            warm=root/f'warm-{name}-{case}.toml'
            if warm.exists():data[name]['long_warmup_median_ms']=record(warm)['median_ms']
        result['storage_diagnostics'][case]=data
    return result


if __name__=='__main__':
    if len(sys.argv)!=3:raise SystemExit(__doc__)
    Path(sys.argv[2]).write_text(json.dumps(summarize(Path(sys.argv[1])),indent=2)+'\n')
