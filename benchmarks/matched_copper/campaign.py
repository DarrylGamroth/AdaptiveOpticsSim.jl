"""Run isolated, serialized validation or repeated service-time measurements."""
import argparse
import json
import os
from pathlib import Path
import random
import resource
import subprocess
import sys
import time
import tomllib

HERE=Path(__file__).resolve().parent
CONFIGS={
    'aos-cpu1':('aos','cpu',1), 'aos-cpu8':('aos','cpu',8),
    'oopao-cpu1':('oopao','cpu',1), 'oopao-cpu8':('oopao','cpu',8),
    'specula-cpu1':('specula','cpu',1),
    'aos-cuda':('aos','cuda',1), 'oopao-cuda':('oopao','cuda',1),
    'specula-cuda':('specula','cuda',1),
}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--mode',choices=('preflight','validate','benchmark'),required=True)
    p.add_argument('--configs',nargs='+',choices=CONFIGS,default=list(CONFIGS))
    p.add_argument('--rounds',type=int)
    args=p.parse_args();root=args.root.resolve()
    protocol=tomllib.loads((HERE/'protocol.toml').read_text())['performance']
    rounds=args.rounds or (protocol['fresh_rounds'] if args.mode=='benchmark' else 1)
    env=os.environ.copy();env.update(JULIA_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONHASHSEED='823')
    order=list(args.configs);random.Random(protocol['order_seed']).shuffle(order)
    journal=root/(args.mode+'-campaign.jsonl')
    for r in range(rounds):
        for key in order[r%len(order):]+order[:r%len(order)]:
            engine,backend,workers=CONFIGS[key]
            output=root/(args.mode+'-'+key+'-r'+str(r+1))
            for candidate in (output,output.with_suffix('.log'),output.with_suffix('.resources.txt')):
                if candidate.exists():raise ValueError('Refusing to overwrite '+str(candidate))
            affinity=(protocol['gpu_affinity'] if backend=='cuda' else
                      protocol['cpu_parallel_affinity'] if workers==8 else protocol['cpu_single_affinity'])
            if engine=='aos':
                cmd=['julia','--project='+str(root/'julia-env'),str(HERE/'run_aos.jl'),
                     str(root/'fixture'),str(output),backend,str(workers),args.mode,str(r+1)]
            else:
                cmd=[str(root/'venv/bin/python'),str(HERE/'run_python.py'),'--engine',engine,
                     '--backend',backend,'--workers',str(workers),'--sources',str(root/'sources'),
                     '--fixture',str(root/'fixture'),'--output',str(output),'--mode',args.mode,'--round',str(r+1)]
                if engine=='specula' and backend=='cuda':cmd.append('--capture')
            cmd=['taskset','-c',affinity]+cmd
            event=dict(config=key,round=r+1,command=cmd,environment={k:env[k] for k in
                ('JULIA_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','PYTHONHASHSEED')},
                started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                loadavg_before=Path('/proc/loadavg').read_text().strip())
            print('Starting',args.mode,key,'round',r+1,flush=True)
            before=resource.getrusage(resource.RUSAGE_CHILDREN);process_start=time.monotonic()
            with output.with_suffix('.log').open('x') as log:
                result=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
            after=resource.getrusage(resource.RUSAGE_CHILDREN)
            event['whole_process_resource_usage']=dict(elapsed_seconds=time.monotonic()-process_start,
                user_cpu_seconds=after.ru_utime-before.ru_utime,system_cpu_seconds=after.ru_stime-before.ru_stime,
                minor_faults=after.ru_minflt-before.ru_minflt,major_faults=after.ru_majflt-before.ru_majflt,
                voluntary_context_switches=after.ru_nvcsw-before.ru_nvcsw,
                involuntary_context_switches=after.ru_nivcsw-before.ru_nivcsw)
            event.update(loadavg_after=Path('/proc/loadavg').read_text().strip(),returncode=result.returncode,finished_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
            with journal.open('a') as f:f.write(json.dumps(event)+'\n')
            if result.returncode:raise RuntimeError('Failed '+key+'; inspect '+str(output.with_suffix('.log')))
    print('Campaign complete',flush=True)

if __name__=='__main__':main()
