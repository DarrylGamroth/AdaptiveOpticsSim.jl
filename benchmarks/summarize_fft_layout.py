"""Hash-check and compare fixed-sampling layout sweeps against native CUDA rates.
Usage: python3 summarize_fft_layout.py EVIDENCE OUTPUT.json
"""
import hashlib
import json
from pathlib import Path
import sys
import tomllib
import numpy as np


def summarize(root):
    result = {'scope': 'Unpaced completed captured physical rate nodes; logical traffic is not measured DRAM traffic', 'cases': {}}
    for case in ('gpi_fine', 'copper_fine'):
        native = tomllib.loads((root / f'native-{case}.toml').read_text())
        first = native['trials'][0]
        reference = np.fromfile(root / first['rate_file'], dtype='<f4').astype(np.float64)
        case_result = {}
        for layout in ('native', 'axis', 'transpose', 'tiled-transpose', 'tiled-pair'):
            path = root / f'{layout}-{case}.toml'
            if not path.exists():
                continue
            meta = tomllib.loads(path.read_text())
            assert meta.get('complete'), path
            batches = {}
            for trial in meta['trials']:
                rate_path = root / trial['rate_file']
                assert hashlib.sha256(rate_path.read_bytes()).hexdigest() == trial['rate_sha256']
                rate = np.fromfile(rate_path, dtype='<f4').astype(np.float64)
                relative = float(np.linalg.norm(rate-reference) / np.linalg.norm(reference))
                flux = float(abs(rate.sum()/reference.sum()-1))
                assert relative < 2e-6 and flux < 2e-6
                assert all(x == 0 for x in trial['host_bytes'])
                batch = str(trial['batch'])
                entry = batches.setdefault(batch, {'round_medians_ms': [], 'maximum_native_relative_l2': 0.0,
                    'maximum_native_flux_relative_error': 0.0, 'field_stack_bytes': trial['field_stack_bytes'],
                    'derived_extra_array_bytes': trial['field_stack_bytes'] * {'native':0,'axis':0,'transpose':2,'tiled-transpose':2,'tiled-pair':1}[layout]})
                entry['round_medians_ms'].append(trial['median_ms'])
                entry['maximum_native_relative_l2'] = max(entry['maximum_native_relative_l2'], relative)
                entry['maximum_native_flux_relative_error'] = max(entry['maximum_native_flux_relative_error'], flux)
                if 'extra_fft_array_bytes' in trial:
                    assert trial['extra_fft_array_bytes'] == entry['derived_extra_array_bytes']
            for entry in batches.values():
                entry['median_of_round_medians_ms'] = float(np.median(entry['round_medians_ms']))
            case_result[layout] = batches
        result['cases'][case] = case_result
    return result


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    Path(sys.argv[2]).write_text(json.dumps(summarize(Path(sys.argv[1])), indent=2)+'\n')
