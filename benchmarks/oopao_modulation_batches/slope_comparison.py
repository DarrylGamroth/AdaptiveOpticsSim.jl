import os,sys,json,hashlib,subprocess
from pathlib import Path
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MPLBACKEND']='Agg'
import numpy as np
root=Path(sys.argv[1] if len(sys.argv)>1 else os.environ.get('AOS_ORACLE_OUTPUT',Path(__file__).parent))
checkout=Path(os.environ.get('OOPAO_PATH',root.parent/'OOPAO'))
sys.path.insert(0,str(checkout))
from OOPAO.Telescope import Telescope
from OOPAO.Source import Source
from OOPAO.Pyramid import Pyramid
p=root/'n120_p30_m3'
tel=Telescope(resolution=120,diameter=8.,samplingTime=.001,centralObstruction=0.,display_optical_path=False)
src=Source(optBand='I',magnitude=0.);src**tel
w=Pyramid(nSubap=30,telescope=tel,modulation=3.,lightRatio=0.,postProcessing='slopesMaps',psfCentering=True,n_pix_separation=2,n_pix_edge=1,binning=1,nTheta_user_defined=32)
def linear_parts(frame):
 q=[w.grabQuadrant(i,cameraFrame=frame)*w.validI4Q for i in range(1,5)]
 a,b,c,d=q
 return np.concatenate((a-b+d-c,a-d+b-c)),np.mean((a+b+c+d)[w.validI4Q])
def signal(frame):
 # The external routine intentionally reads camera storage, ignoring its frame argument.
 w.cam.frame=frame.copy()
 return w.signalProcessing()[1].copy()
def tangent(flat,derivative):
 s,n=linear_parts(flat);ds,dn=linear_parts(derivative)
 return ((ds/n-s*dn/n**2)*w.slopesUnits)[w.validSignal==1]
def load(name):return np.loadtxt(p/(name+'.txt'))
def metrics(a,b):
 gain=float(np.dot(a,b)/np.dot(b,b))
 return dict(relative_l2=float(np.linalg.norm(a-b)/np.linalg.norm(b)),gain=gain,cosine=float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b))),orthogonal_relative_l2=float(np.linalg.norm(a-gain*b)/np.linalg.norm(b)))
f=load('flat_frame');r=load('ramp_frame');d=load('tip_frame_derivative')
ref_tip=load('tip_signal_derivative');ref_ramp=load('ramp_signal')-load('flat_signal')
report={'validation':{'saved_ramp_signal':metrics(signal(r)-signal(f),ref_ramp),'analytic_tip_vs_saved_central_difference':metrics(tangent(f,d),ref_tip)}}
for strategy in ('pupil_tilt','shifted_mask'):
 for batch in ('baseline','selected'):
  prefix=f'cuda_f32_{strategy}_{batch}'
  ff=load(prefix+'_flat_frame');rr=load(prefix+'_ramp_frame');dd=load(prefix+'_tip_frame_derivative')
  report[prefix]={'ramp_slopes':metrics(signal(rr)-signal(ff),ref_ramp),'tip_slopes':metrics(tangent(ff,dd),ref_tip),'tip_image':metrics(dd.ravel(),d.ravel())}
files=['precision_oopao.npy','OOPAO/__init__.py','OOPAO/Pyramid.py','OOPAO/Telescope.py','OOPAO/Source.py','OOPAO/Detector.py']
report['oopao_provenance']={'head':subprocess.check_output(['git','-C',str(checkout),'rev-parse','HEAD'],text=True).strip(),'status':subprocess.check_output(['git','-C',str(checkout),'status','--short'],text=True),'sha256':{f:hashlib.sha256((checkout/f).read_bytes()).hexdigest() for f in files},'precision_value':int(np.load(checkout/'precision_oopao.npy')),'precision_write':'OOPAO/__init__.py unconditionally writes 64 on import; no restoration attempted'}
assert report['validation']['saved_ramp_signal']['relative_l2'] < 1e-12
assert report['validation']['analytic_tip_vs_saved_central_difference']['relative_l2'] < 1e-6
(root/'slope-comparison.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
