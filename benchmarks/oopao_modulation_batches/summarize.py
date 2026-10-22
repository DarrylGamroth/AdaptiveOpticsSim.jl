from pathlib import Path
import numpy as np,json,subprocess,os,sys
root=Path(sys.argv[1] if len(sys.argv)>1 else os.environ.get('AOS_ORACLE_OUTPUT',Path(__file__).parent))
prefix=os.environ.get('AOS_ORACLE_PREFIX','aos')
report=json.loads((root/'summary.json').read_text())
for name in ('n24_p4_m1','n32_p8_m3','n120_p30_m3'):
 if name not in os.environ.get('AOS_ORACLE_CASES','n24_p4_m1,n32_p8_m3,n120_p30_m3').split(','): continue
 p=root/name;flat_ref=np.loadtxt(p/'flat_frame.txt');flat=np.loadtxt(p/(prefix+'_flat_frame.txt'))
 metrics={}
 for label in ('flat_frame','ramp_frame','tip_frame_derivative'):
  ref=np.loadtxt(p/(label+'.txt'));a=np.loadtxt(p/(prefix+'_'+label+'.txt'))
  metrics[label]={'relative_l2':float(np.linalg.norm(a-ref)/np.linalg.norm(ref)), 'flux_normalized_relative_l2':float(np.linalg.norm(a/flat.sum()-ref/flat_ref.sum())/np.linalg.norm(ref/flat_ref.sum()))}
 metrics['flat_photon_flux_relative_error']=float(abs(flat.sum()-flat_ref.sum())/flat_ref.sum())
 metrics['aos_flat_photons']=float(flat.sum());report[name][prefix+'_comparison']=metrics
if os.environ.get('AOS_CHECKOUT'):
 report['aos_commit']=subprocess.check_output(['git','-C',os.environ['AOS_CHECKOUT'],'rev-parse','HEAD'],text=True).strip()
report['aos_worktree']=os.environ.get('AOS_CHECKOUT','unspecified; record Julia environment manifest')
report['comparison']='OOPAO physical pupil tilt versus AOS '+prefix+' '+os.environ.get('AOS_ORACLE_STRATEGY','pupil_tilt')+'; shifted-mask comparison measures an existing approximation'
(root/('comparison.json' if prefix=='aos' else 'comparison_'+prefix+'.json')).write_text(json.dumps(report,indent=2))
