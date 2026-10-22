"""Native OOPAO/SPECULA optical execution under the fixed common prescription."""
from pathlib import Path
import argparse
import hashlib
import io
import json
import os
import platform
import subprocess
import sys
import time
import tomllib
from contextlib import redirect_stdout

import numpy as np
from common import load_fixture,save_array,write_json


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--engine',choices=('oopao','specula'),required=True)
    p.add_argument('--backend',choices=('cpu','cuda'),required=True)
    p.add_argument('--workers',type=int,default=1)
    p.add_argument('--capture',action='store_true')
    p.add_argument('--native',action='store_true')
    p.add_argument('--sources',type=Path,required=True)
    p.add_argument('--fixture',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mode',choices=('preflight','validate','benchmark'),default='preflight')
    p.add_argument('--round',type=int,default=0)
    args=p.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('Output must be empty')
    args.output.mkdir(parents=True,exist_ok=True)
    src=args.sources/args.engine
    lock=json.loads((Path(__file__).parent/'sources.json').read_text())
    expected=next(s['revision'] for s in lock['sources'] if s['name'].lower()==args.engine)
    head=subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip()
    if head!=expected:
        raise ValueError('External source revision differs from sources.json')
    dirty=subprocess.check_output(['git','-C',str(src),'status','--porcelain'],text=True)
    if dirty:
        raise ValueError('External source checkout is dirty: '+dirty)
    sys.path.insert(0,str(src))
    fixture,arrays=load_fixture(args.fixture)
    contract=fixture['contract'];probes=arrays['probes']
    gpu=args.backend=='cuda'
    if args.engine=='oopao':
        if args.native or args.capture:
            raise ValueError('OOPAO uses its native batched propagation; capture is not supplied by this adapter')
        os.environ.update(OOPAO_BACKEND=args.backend,OOPAO_PRECISION='32',OOPAO_GPU_RESIDENT='1' if gpu else '0')
        from oopao_adapter import OopaoRunner
        make=lambda:OopaoRunner(contract,arrays['amplitude'],probes[:,:,0],backend=args.backend,cpu_workers=args.workers)
    else:
        import specula
        specula.init(0 if gpu else -1,precision=1)
        from specula_adapter import CopperSpeculaRunner
        make=lambda:CopperSpeculaRunner(contract,arrays['amplitude'],probes[:,:,0],
                                        target_device_idx=0 if gpu else -1,precision=1,
                                        matched=not args.native,capture=args.capture)
    if gpu:
        import cupy as cp
        wait=lambda:cp.cuda.runtime.deviceSynchronize()
        host=cp.asnumpy
    else:
        wait=lambda:None
        host=np.asarray
    start=time.perf_counter_ns();runner=make();wait();preparation=(time.perf_counter_ns()-start)/1e9
    def stage(opd):
        runner.upload_opd(opd)
        wait()
    def step():
        value=runner.step()
        wait()
        return value
    rates=np.empty((64,64,probes.shape[2]),dtype=np.float32)
    for i in range(probes.shape[2]):
        stage(probes[:,:,i]);value=step()
        if tuple(value.shape)!=(64,64) or value.dtype!=np.float32:
            raise ValueError('Output shape/precision mismatch')
        if gpu and not isinstance(value,cp.ndarray):
            raise ValueError('CUDA output is not resident')
        rates[:,:,i]=host(value)
    if not np.all(np.isfinite(rates)) or np.min(rates)<0:
        raise ValueError('Invalid photon rates')
    if np.array_equal(rates[:,:,0],rates[:,:,7]):
        raise ValueError('Stale captured input')
    protocol=tomllib.loads((Path(__file__).parent/'protocol.toml').read_text())
    files=['run_python.py','common.py',args.engine+'_adapter.py','protocol.toml','sources.json']
    record=dict(engine=args.engine,backend=args.backend,cpu_workers=args.workers,mode=args.mode,
                round=args.round,complete=False,prescription='native-diagnostic' if args.native else 'common',
                external_revision=head,external_status_before=dirty,python=platform.python_version(),
                numpy=np.__version__,preparation_seconds=preparation,metadata=runner.metadata(),
                source_sha256={f:hashlib.sha256((Path(__file__).parent/f).read_bytes()).hexdigest() for f in files},
                fixture_sha256=hashlib.sha256((args.fixture/'fixture.toml').read_bytes()).hexdigest(),
                rates=save_array(args.output,'rates',rates),boundary=protocol['scope'],arrival=protocol['arrival'])
    if gpu:
        config=io.StringIO()
        with redirect_stdout(config):cp.show_config()
        record.update(cupy=cp.__version__,cuda_configuration=config.getvalue(),
                      cufft_version=cp.cuda.cufft.getVersion())
    write_json(args.output/'record.json',record)
    if args.mode=='validate':
        basis=arrays['calibration_basis'];poke=np.float32(contract['calibration_poke_m'])
        positive=np.empty((64,64,277),dtype=np.float32);negative=np.empty_like(positive)
        scratch=np.empty((240,240),dtype=np.float32)
        for k in range(277):
            np.multiply(basis[:,:,k],poke,out=scratch)
            stage(scratch);positive[:,:,k]=host(step())
            np.multiply(basis[:,:,k],-poke,out=scratch)
            stage(scratch);negative[:,:,k]=host(step())
            if (k+1)%32==0:print('calibration',k+1,'/277',flush=True)
        record.update(calibration_positive=save_array(args.output,'calibration_positive',positive),
                      calibration_negative=save_array(args.output,'calibration_negative',negative))
        write_json(args.output/'record.json',record)
    elif args.mode=='benchmark':
        perf=protocol['performance'];warmups=perf['warmups'];count=perf['samples']
        for k in range(warmups):
            stage(probes[:,:,k%probes.shape[2]]);step()
        memory_before=cp.get_default_memory_pool().used_bytes() if gpu else None
        samples=[]
        for k in range(count):
            stage(probes[:,:,k%probes.shape[2]])
            start=time.perf_counter_ns();value=step();samples.append((time.perf_counter_ns()-start)/1e6)
        if not np.array_equal(host(value),rates[:,:,(count-1)%probes.shape[2]]):
            raise ValueError('Repeated resident OPD replay changed output')
        memory_after=cp.get_default_memory_pool().used_bytes() if gpu else None
        record.update(service=dict(samples_ms=samples,median_ms=float(np.median(samples)),
                                    minimum_ms=min(samples),maximum_ms=max(samples),
                                    p95_ms=float(np.quantile(samples,.95)),warmups=warmups,samples=count),
                      observed_default_pool_before=memory_before,observed_default_pool_after=memory_after)
    record['complete']=True
    record['external_status_after']=subprocess.check_output(['git','-C',str(src),'status','--porcelain'],text=True)
    write_json(args.output/'record.json',record)


if __name__=='__main__':
    main()
