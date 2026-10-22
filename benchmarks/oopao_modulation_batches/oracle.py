import sys,os,json,subprocess
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MPLBACKEND']='Agg'
import numpy as np
from pathlib import Path
checkout=Path(os.environ.get('OOPAO_PATH',Path(__file__).resolve().parent.parent/'OOPAO')).resolve()
sys.path.insert(0,str(checkout))
from OOPAO.Telescope import Telescope
from OOPAO.Source import Source
from OOPAO.Pyramid import Pyramid
from pathlib import Path
root=Path(sys.argv[1] if len(sys.argv)>1 else os.environ.get('AOS_ORACLE_OUTPUT',Path(__file__).parent))
root.mkdir(parents=True,exist_ok=True)
summary=json.loads((root/'summary.json').read_text()) if (root/'summary.json').exists() else {}
for resolution,pup,mod,points,sep,edge,binning in [(24,4,1.,8,4,0,2),(32,8,3.,20,8,0,1),(120,30,3.,32,2,1,1)]:
 name=f'n{resolution}_p{pup}_m{int(mod)}'
 if name not in os.environ.get('AOS_ORACLE_CASES','n24_p4_m1,n32_p8_m3,n120_p30_m3').split(','): continue
 out=root/name;out.mkdir(exist_ok=True)
 tel=Telescope(resolution=resolution,diameter=8.,samplingTime=.001,centralObstruction=0.,display_optical_path=False)
 src=Source(optBand='I',magnitude=0.);src**tel
 wfs=Pyramid(nSubap=pup,telescope=tel,modulation=mod,lightRatio=0.,postProcessing='slopesMaps',psfCentering=True,n_pix_separation=sep,n_pix_edge=edge,binning=binning,nTheta_user_defined=points)
 def save(label,opd):
  tel.OPD_no_pupil=opd;tel*wfs
  np.savetxt(out/(label+'_frame.txt'),wfs.cam.frame)
  np.savetxt(out/(label+'_raw.txt'),wfs.raw_data)
  np.savetxt(out/(label+'_signal.txt'),wfs.signal)
  return np.array(wfs.cam.frame),np.array(wfs.signal)
 zero=np.zeros((resolution,resolution)); f,s=save('flat',zero)
 x=np.arange(resolution)[:,None];y=np.arange(resolution)[None,:]
 ramp=5e-9*x-2e-9*y
 fr,sr=save('ramp',ramp)
 basis=np.broadcast_to(np.linspace(-1.,1.,resolution,endpoint=False)[:,None],zero.shape)*tel.pupil
 plus,sp=save('tip_plus',1e-9*basis);minus,sm=save('tip_minus',-1e-9*basis)
 np.savetxt(out/'tip_frame_derivative.txt',(plus-minus)/2e-9);np.savetxt(out/'tip_signal_derivative.txt',(sp-sm)/2e-9)
 summary[name]=dict(resolution=resolution,pupil_samples=pup,modulation=mod,modulation_points=wfs.nTheta,n_pix_separation=sep,n_pix_edge=edge,binning=binning,flat_flux=float(f.sum()),ramp_flux=float(fr.sum()),source_flux_per_frame=float(src.fluxMap.sum()),wavelength=float(src.wavelength),flat_frame_shape=list(f.shape),raw_shape=list(wfs.raw_data.shape),flat_signal_norm=float(np.linalg.norm(s)),tip_derivative_norm=float(np.linalg.norm((sp-sm)/2e-9)))
summary['oopao_commit']=subprocess.check_output(['git','-C',str(checkout),'rev-parse','HEAD'],text=True).strip()
(root/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
