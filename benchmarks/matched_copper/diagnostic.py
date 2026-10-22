"""CPU-only cold-factor diagnostic; never produces accepted benchmark records.

Exports realized AOS/OOPAO centering and prism arrays. On flat and 20 nm x-tilt
only, substitutes OOPAO cold factors into a separate AOS instance one at a time
and together. This is an intervention to localize numerical differences, not a
physical correction, production change, tolerance adjustment, or timing run.
"""
import os
os.environ.update(OOPAO_BACKEND='cpu', OOPAO_PRECISION='32', OOPAO_GPU_RESIDENT='0',
                  OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import numpy as np
from common import load_fixture, incident_divisor, slopes, write_json

JULIA_DIAGNOSTIC = r'''
include(ARGS[1])
function diagnostic(fixture, output)
    BLAS.set_num_threads(1)
    Backends.set_fft_provider_threads!(1)
    meta=TOML.parsefile(joinpath(fixture,"fixture.toml"))
    graph=make_graph(meta["contract"],false)
    owner=AlgorithmGraphs.prepared_graph_node(graph,:pwfs).prepared
    Array(owner.input.amplitude)==read_array(fixture,meta["amplitude"]) || error("pupil mismatch")
    workspace=owner.workspace
    original_phasor=copy(workspace.phasor)
    original_mask=copy(workspace.pyramid_mask)
    for (name,array) in (("aos_phasor",original_phasor),("aos_mask",original_mask))
        open(io->write(io,vec(array)),joinpath(output,name*".c64"),"w")
    end
    foreign_phasor=reshape(copy(reinterpret(ComplexF32,read(joinpath(output,"oopao_phasor.c64")))),512,512)
    foreign_mask=reshape(copy(reinterpret(ComplexF32,read(joinpath(output,"oopao_mask.c64")))),512,512)
    probes=read_array(fixture,meta["probes"])
    input=graph_input(graph,:opd)
    rate=graph_output(graph,:rate)
    for (name,phasor,mask) in (("baseline",original_phasor,original_mask),
        ("phasor_only",foreign_phasor,original_mask),
        ("mask_only",original_phasor,foreign_mask),
        ("both",foreign_phasor,foreign_mask))
        copyto!(workspace.phasor,phasor)
        copyto!(workspace.pyramid_mask,mask)
        result=Array{Float32}(undef,64,64,2)
        for (slot,index) in enumerate((1,3))
            copyto!(input,view(probes,:,:,index))
            step_graph!(graph)
            result[:,:,slot].=rate
        end
        open(io->write(io,vec(result)),joinpath(output,"aos_"*name*".f32"),"w")
    end
end
diagnostic(ARGS[2],ARGS[3])
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources',type=Path,required=True)
    parser.add_argument('--fixture',type=Path,required=True)
    parser.add_argument('--julia-project',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('Diagnostic output must be empty')
    args.output.mkdir(parents=True,exist_ok=True)
    here=Path(__file__).resolve().parent
    lock=json.loads((here/'sources.json').read_text())
    oopao=next(item for item in lock['sources'] if item['name']=='OOPAO')
    source=args.sources/oopao['path']
    revision=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    status=subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True)
    if revision!=oopao['revision'] or status:
        raise ValueError('OOPAO source must match the clean pinned revision')
    sys.path.insert(0,str(source))
    from oopao_adapter import OopaoRunner
    fixture,arrays=load_fixture(args.fixture)
    runner=OopaoRunner(fixture['contract'],arrays['amplitude'],arrays['probes'][:,:,0],cpu_workers=1)
    for name,array in (('phasor',runner.sensor.phasor),('mask',runner.sensor.mask)):
        values=np.asarray(array,dtype='<c8')
        if values.shape!=(512,512):raise ValueError('Unexpected cold array shape')
        (args.output/('oopao_'+name+'.c64')).write_bytes(values.tobytes(order='F'))
    external=np.empty((64,64,2),dtype=np.float32)
    for slot,index in enumerate((0,2)):
        runner.upload_opd(arrays['probes'][:,:,index])
        external[:,:,slot]=runner.step()
    (args.output/'oopao_rates.f32').write_bytes(external.tobytes(order='F'))
    subprocess.run(['julia','--threads=1','--project='+str(args.julia_project),'-e',JULIA_DIAGNOSTIC,
                    str(here/'run_aos.jl'),str(args.fixture),str(args.output)],check=True)
    cold={}
    for name in ('phasor','mask'):
        reference=np.fromfile(args.output/('oopao_'+name+'.c64'),dtype='<c8').astype(np.complex128)
        candidate=np.fromfile(args.output/('aos_'+name+'.c64'),dtype='<c8').astype(np.complex128)
        error=candidate-reference
        cold[name]=dict(maximum_absolute=float(np.max(np.abs(error))),
                        relative_l2=float(np.linalg.norm(error)/np.linalg.norm(reference)),
                        maximum_absolute_phase_rad=float(np.max(np.abs(np.angle(candidate*np.conj(reference))))))
    divisor=incident_divisor(fixture['contract'],arrays['amplitude'])
    variants={}
    for variant in ('baseline','phasor_only','mask_only','both'):
        values=np.fromfile(args.output/('aos_'+variant+'.f32'),dtype='<f4').reshape((64,64,2),order='F')
        comparisons=[]
        for slot,name in enumerate(('flat','tilt_x_20nm')):
            reference=external[:,:,slot].astype(float);candidate=values[:,:,slot].astype(float)
            rs=slopes(reference,divisor);cs=slopes(candidate,divisor)
            comparisons.append(dict(probe=name,rate_relative_l2=float(np.linalg.norm(candidate-reference)/np.linalg.norm(reference)),
                slope_absolute_max=float(np.max(np.abs(cs-rs))),slope_error_l2=float(np.linalg.norm(cs-rs)),
                slope_relative_l2=float(np.linalg.norm(cs-rs)/np.linalg.norm(rs)),
                reference_slope_l2=float(np.linalg.norm(rs))))
        variants[variant]=comparisons
    files={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in args.output.iterdir() if path.is_file()}
    result=dict(scope='CPU-only two-probe cold-factor intervention; no qualification or timing record',
                oopao_revision=revision,probes=['flat','tilt_x_20nm'],cold=cold,variants=variants,
                files_sha256=files,source_sha256={name:hashlib.sha256((here/name).read_bytes()).hexdigest()
                for name in ('diagnostic.py','run_aos.jl','oopao_adapter.py','common.py','protocol.toml')},
                fixture_sha256=hashlib.sha256((args.fixture/'fixture.toml').read_bytes()).hexdigest())
    write_json(args.output/'diagnostic.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
