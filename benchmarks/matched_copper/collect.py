"""Reduce preserved numerical checks and three fresh performance rounds."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tomllib
import numpy as np
from common import write_json
from evaluate import digest,load_run,compare_arrays
from common import load_fixture,incident_divisor,support_indices
from campaign import CONFIGS

HERE=Path(__file__).resolve().parent
VALIDATION={
 'aos-cpu1':'aos-cpu-validate','aos-cpu8':'validate-aos-cpu8-r1',
 'oopao-cpu1':'oopao-cpu-validate','oopao-cpu8':'validate-oopao-cpu8-r1',
 'specula-cpu1':'specula-cpu-validate','aos-cuda':'aos-cuda-validate',
 'oopao-cuda':'oopao-cuda-validate','specula-cuda':'specula-cuda-validate',
}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
 p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.root.resolve()
 if a.output.exists() and any(a.output.iterdir()):raise ValueError('Output must be empty')
 a.output.mkdir(parents=True,exist_ok=True)
 for key in VALIDATION:
  standard='validate-'+key+'-r1'
  if (root/standard).exists():VALIDATION[key]=standard
 fixture,fixture_arrays=load_fixture(root/'fixture')
 divisor=incident_divisor(fixture['contract'],fixture_arrays['amplitude'])
 incident=divisor*len(support_indices()[0])
 fixture_hash=digest(root/'fixture/fixture.toml');protocol=tomllib.loads((HERE/'protocol.toml').read_text())
 protocol_hash=digest(HERE/'protocol.toml');checks={}
 pairs=[('cpu-oopao','aos-cpu1','oopao-cpu1',False),('cpu-specula','aos-cpu1','specula-cpu1',False),
        ('cuda-oopao','aos-cuda','oopao-cuda',False),('cuda-specula','aos-cuda','specula-cuda',False),
        ('aos-cpu-cuda','aos-cpu1','aos-cuda',True),('oopao-cpu-cuda','oopao-cpu1','oopao-cuda',True),
        ('specula-cpu-cuda','specula-cpu1','specula-cuda',True),
        ('aos-cpu-threading','aos-cpu1','aos-cpu8',False),('oopao-cpu-threading','oopao-cpu1','oopao-cpu8',False)]
 for name,ref,can,cross in pairs:
  out=a.output/name
  cmd=[sys.executable,str(HERE/'evaluate.py'),'--reference',str(root/VALIDATION[ref]),
       '--candidate',str(root/VALIDATION[can]),'--fixture',str(root/'fixture'),'--output',str(out)]
  if cross:cmd.append('--cpu-cuda')
  result=subprocess.run(cmd,capture_output=True,text=True)
  (a.output/(name+'.log')).write_text(result.stdout+result.stderr)
  if not (out/'summary.json').exists():raise RuntimeError(name+' evidence could not be evaluated')
  d=json.loads((out/'summary.json').read_text());cal=d['calibration']
  checks[name]=dict(passed=d['passed'],full_validation=d['full_validation'],timing_eligible=d['timing_eligible'],
      maximum_rate_relative_l2=max(v['rate_relative_l2'] for v in d['probes']),
      maximum_flux_relative_error=max(v['flux_relative_error'] for v in d['probes']),
      maximum_raw_slope_absolute_error=max(v['slope_absolute_error'] for v in d['probes']),
      maximum_excited_slope_relative_l2=max(v['slope_relative_l2'] for v in d['probes'][2:]),
      calibration_relative_l2=cal['relative_l2'],maximum_column_relative_l2=cal['maximum_column_relative_l2'],
      failed_checks=[v for v in d['checks'] if not v['passed']])
 performance={}
 for key in CONFIGS:
  validation,va,_=load_run(root/VALIDATION[key],fixture_hash,protocol_hash,fixture['aos_revision']);rounds=[]
  expected=CONFIGS[key]
  if tuple(validation[k] for k in ('engine','backend','cpu_workers'))!=expected:
   raise ValueError('Wrong validation configuration: '+key)
  for r in range(1,protocol['performance']['fresh_rounds']+1):
   rec,arrays,_=load_run(root/('benchmark-'+key+'-r'+str(r)),fixture_hash,protocol_hash,fixture['aos_revision'])
   if tuple(rec[k] for k in ('engine','backend','cpu_workers'))!=expected:
    raise ValueError('Wrong timing configuration: '+key)
   if rec['mode']!='benchmark' or rec['round']!=r:raise ValueError('Wrong performance round')
   if rec['source_sha256']!=validation['source_sha256']:raise ValueError('Different qualification implementation')
   replay,_=compare_arrays({'rates':va['rates']},{'rates':arrays['rates']},
       fixture['probe_names'],divisor,incident,fixture['contract']['calibration_poke_m'],
       protocol['validation'],cpu_cuda=True)
   if not replay['passed']:raise ValueError('Validation/performance numerical parity failed: '+key)
   if rec.get('actual_batch')!=validation.get('actual_batch') or rec.get('metadata')!=validation.get('metadata'):
    raise ValueError('Different validation/performance configuration: '+key)
   series=rec['service']['samples_ms']
   if len(series)!=protocol['performance']['samples'] or not np.all(np.isfinite(series)) or min(series)<=0:
    raise ValueError('Invalid timing series')
   reduced=dict(samples_ms=series,median_ms=float(np.median(series)),minimum_ms=min(series),
       maximum_ms=max(series),p95_ms=float(np.quantile(series,.95)),
       warmups=rec['service']['warmups'],samples=len(series))
   if reduced['warmups']!=protocol['performance']['warmups']:raise ValueError('Wrong warmup count')
   for metric in ('median_ms','minimum_ms','maximum_ms','p95_ms'):
    if not np.isclose(reduced[metric],rec['service'][metric],rtol=1e-12,atol=1e-12):
     raise ValueError('Incorrect saved statistic: '+metric)
   rounds.append(dict(round=r,**reduced,validation_probe_parity=dict(passed=replay['passed'],
       maximum_rate_relative_l2=max(v['rate_relative_l2'] for v in replay['probes']),
       maximum_raw_slope_absolute_error=max(v['slope_absolute_error'] for v in replay['probes']),
       rate_limit=protocol['validation']['cpu_cuda_rate_relative_l2']),preparation_seconds=rec['preparation_seconds'],
       host_allocation_bytes=rec.get('host_allocation_bytes'),
       observed_cuda_pool_bytes=rec.get('observed_cuda_pool_bytes'),
       observed_default_pool_before=rec.get('observed_default_pool_before'),
       observed_default_pool_after=rec.get('observed_default_pool_after')))
  medians=[r['median_ms'] for r in rounds]
  performance[key]=dict(median_of_round_medians_ms=float(np.median(medians)),
      round_median_range_ms=[min(medians),max(medians)],rounds=rounds,
      metadata=validation.get('metadata'),actual_batch=validation.get('actual_batch'),
      interpretation=('Descriptive shared-host optical service time; numerical qualification '+
                      ('passed.' if all(c['passed'] for c in checks.values()) else 'has failed gates; see numerical results.')))
 write_json(a.output/'summary.json',dict(numerical=checks,performance=performance,
     source_revision=tomllib.loads((root/'fixture/fixture.toml').read_text())['aos_revision'],
     fixture_sha256=fixture_hash,protocol_sha256=protocol_hash,
     collector_sha256=digest(__file__),scope=protocol['scope'],arrival=protocol['arrival']))
 print('Saved',a.output/'summary.json')

if __name__=='__main__':main()
