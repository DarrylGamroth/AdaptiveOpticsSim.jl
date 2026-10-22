"""Compare fine FFT implementation changes against fixed-sampling baselines.
Usage: python3 compare_fft_validation.py PREVIOUS_REUSE_EVIDENCE FFT_EVIDENCE OUTPUT.json
"""
import json
from pathlib import Path
import sys
import tomllib
from compare_pupil_reuse_validation import compare_copper, compare_gpi


def compare(previous, current):
    old_dir, new_dir = current/'copper-before', current/'copper-after'
    for directory, expected_batch in ((old_dir,2),(new_dir,16)):
        record = tomllib.loads((directory/'validation.toml').read_text())['cases']['n1440-m256']
        assert record['resolution'] == 1440 and record['points'] == 256
        assert record['modulation_batch']['batch_size'] == expected_batch
    preflight = tomllib.loads((current/'gpi-after.preflight.toml').read_text())
    assert preflight['fine']['batch_size'] == 16
    assert 'PyramidTransposeFFTWorkspace' in preflight['fine']['fft_workspace_type']
    result = {'scope':'Fixed physical sampling; fine COPPER frozen calibration and complete captured GPI loops',
              'copper':compare_copper(old_dir,new_dir,key='n1440-m256'),
              'gpi':compare_gpi(previous/'gpi-after',current/'gpi-after'), 'passed':True}
    return result


if __name__ == '__main__':
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    Path(sys.argv[3]).write_text(json.dumps(compare(Path(sys.argv[1]),Path(sys.argv[2])),indent=2)+'\n')
