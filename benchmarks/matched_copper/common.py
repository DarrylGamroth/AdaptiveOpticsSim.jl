"""Shared file format and readout, outside all engine timing boundaries."""
from pathlib import Path
import hashlib
import json
import tomllib
import numpy as np


def load_fixture(path):
    root=Path(path)
    record=tomllib.loads((root/'fixture.toml').read_text())
    arrays={key:read_array(root,record[key]) for key in ('amplitude','probes','calibration_basis')}
    return record,arrays


def read_array(root, meta):
    data=(Path(root)/meta['file']).read_bytes()
    if hashlib.sha256(data).hexdigest()!=meta['sha256']:
        raise ValueError('Array digest mismatch: '+meta['file'])
    return np.frombuffer(data,dtype='<f4').reshape(meta['shape'],order='F')


def save_array(root,name,array):
    values=np.asarray(array,dtype='<f4')
    path=Path(root)/(name+'.f32')
    path.write_bytes(values.tobytes(order='F'))
    return dict(file=path.name,shape=list(values.shape),dtype='float32',order='F',
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def support_indices():
    x,y=np.meshgrid(np.arange(30)-14.5,np.arange(30)-14.5,indexing='ij')
    flat=np.flatnonzero((x*x+y*y<=225).ravel(order='F'))
    return flat%30,flat//30


def incident_divisor(contract,amplitude):
    # Physical incident rate / number of supported samples. This is fixed
    # analytically and shared; no external image or response is gain-fitted.
    pitch=float(np.float32(contract['diameter_m'])/np.float32(contract['resolution']))
    incident=float(contract['photon_irradiance'])*pitch*pitch*np.sum(amplitude.astype(float)**2)
    return incident/len(support_indices()[0])


def slopes(rate,divisor):
    i,j=support_indices()
    a,b,c,d=(np.asarray(rate,dtype=float)[col+j,row+i]
             for row,col in ((1,1),(33,1),(33,33),(1,33)))
    return np.concatenate(((a-b+d-c)/divisor,(a-d+b-c)/divisor))


def write_json(path,record):
    path=Path(path)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)
