"""Evaluate raw common-prescription rates and central calibration evidence.

No gain, flat subtraction, image transform, or calibration column rescaling is
performed. Derived slopes and interaction matrices retain Float64 precision.
"""

from pathlib import Path
import argparse
import hashlib
import json
import sys
import tomllib

import numpy as np

from common import incident_divisor, load_fixture, read_array, slopes, write_json


HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_f64(root, name, values):
    path = Path(root) / (name + '.f64')
    values = np.asarray(values, dtype='<f8')
    path.write_bytes(values.tobytes(order='F'))
    return dict(file=path.name, shape=list(values.shape), dtype='float64',
                order='F', sha256=digest(path))


def load_run(root, fixture_hash, protocol_hash, aos_revision):
    root = Path(root)
    paths = [p for p in (root/'record.toml', root/'record.json') if p.exists()]
    if len(paths) != 1:
        raise ValueError(f'{root}: expected exactly one record.toml or record.json')
    path = paths[0]
    record = (tomllib.loads(path.read_text()) if path.suffix == '.toml'
              else json.loads(path.read_text()))
    if not record.get('complete'):
        raise ValueError(f'{root}: run is incomplete')
    if record['fixture_sha256'] != fixture_hash:
        raise ValueError(f'{root}: fixture hash differs')
    if record['source_sha256']['protocol.toml'] != protocol_hash:
        raise ValueError(f'{root}: frozen protocol hash differs')
    engine = record['engine']
    required = ({'run_aos.jl', 'protocol.toml'} if engine == 'aos' else
                {'run_python.py', 'common.py', engine+'_adapter.py', 'protocol.toml', 'sources.json'})
    if not required.issubset(record['source_sha256']):
        raise ValueError(f'{root}: missing harness source hashes')
    for name, expected_hash in record['source_sha256'].items():
        if Path(name).name != name or not (HERE/name).is_file() or digest(HERE/name) != expected_hash:
            raise ValueError(f'{root}: current harness bytes differ for {name}')
    if engine == 'aos':
        expected_revision = aos_revision
    else:
        source_lock = json.loads((HERE/'sources.json').read_text())
        revisions = {source['name'].lower(): source['revision'] for source in source_lock['sources']}
        if engine not in revisions:
            raise ValueError(f'{root}: engine absent from exact source lock: {engine}')
        expected_revision = revisions[engine]
        if record.get('external_status_before') != '' or record.get('external_status_after') != '':
            raise ValueError(f'{root}: external source checkout was not recorded clean')
    if record.get('external_revision', record.get('revision')) != expected_revision:
        raise ValueError(f'{root}: source revision differs from exact source lock/fixture')
    arrays = {'rates': read_array(root, record['rates'])}
    for name in ('calibration_positive', 'calibration_negative'):
        if name in record:
            arrays[name] = read_array(root, record[name])
    if ('calibration_positive' in arrays) != ('calibration_negative' in arrays):
        raise ValueError(f'{root}: incomplete central calibration pair')
    if record['mode'] == 'validate' and 'calibration_positive' not in arrays:
        raise ValueError(f'{root}: validate record lacks calibration arrays')
    for name, array in arrays.items():
        expected = (64, 64, 17 if name == 'rates' else 277)
        meta = record[name]
        if array.shape != expected or meta['dtype'] != 'float32' or meta['order'] != 'F':
            raise ValueError(f'{root}: invalid {name} format; expected Float32 {expected} in F order')
        if not np.isfinite(array).all() or np.min(array) < 0:
            raise ValueError(f'{root}: nonfinite or negative {name}')
    provenance = dict(directory=str(root.resolve()), record_file=path.name,
                      record_sha256=digest(path), engine=record['engine'],
                      backend=record['backend'], mode=record['mode'],
                      prescription=record['prescription'],
                      source_sha256=record['source_sha256'],
                      revision=record.get('external_revision', record.get('revision')),
                      arrays={key: record[key] for key in arrays})
    return record, arrays, provenance


def readout_stack(rates, divisor):
    return np.stack([slopes(rates[:, :, k], divisor)
                     for k in range(rates.shape[2])], axis=1)


def relative_l2(error, reference):
    denominator = float(np.linalg.norm(reference))
    numerator = float(np.linalg.norm(error))
    if denominator == 0:
        return 0.0 if numerator == 0 else None
    return numerator / denominator


def gate(name, value, limit):
    return dict(name=name, value=value, limit=float(limit),
                passed=value is not None and bool(value <= limit))


def compare_arrays(reference, candidate, probe_names, divisor, incident, poke,
                   thresholds, cpu_cuda=False):
    """Pure numeric core; raw slopes are retained, including fixed flat bias."""
    rr, cr = reference['rates'].astype(float), candidate['rates'].astype(float)
    rs, cs = readout_stack(rr, divisor), readout_stack(cr, divisor)
    checks, probes = [], []
    rate_limit = thresholds['cpu_cuda_rate_relative_l2'] if cpu_cuda else thresholds['rate_relative_l2']
    flat = probe_names.index('flat')
    excited = probe_names.index('tilt_x_100nm')
    live_input = {}
    for owner, rates in (('reference', rr), ('candidate', cr)):
        changed = not np.array_equal(rates[:, :, flat], rates[:, :, excited])
        live_input[owner] = dict(probe='tilt_x_100nm', changed=changed,
                                 rate_relative_l2=relative_l2(rates[:, :, excited]-rates[:, :, flat],
                                                             rates[:, :, flat]))
        checks.append(dict(name=owner+':live_input', passed=changed,
                           rule='flat and tilt_x_100nm outputs must differ exactly'))
    for k, name in enumerate(probe_names):
        rate_error = relative_l2(cr[:, :, k]-rr[:, :, k], rr[:, :, k])
        slope_error = cs[:, k]-rs[:, k]
        slope_absolute = float(np.max(np.abs(slope_error)))
        slope_relative = relative_l2(slope_error, rs[:, k])
        rf, cf = float(np.sum(rr[:, :, k])), float(np.sum(cr[:, :, k]))
        flux_error = abs(cf-rf)/abs(rf) if rf != 0 else None
        conservation = dict(reference=abs(rf-incident)/incident,
                            candidate=abs(cf-incident)/incident)
        entry = dict(name=name, rate_relative_l2=rate_error,
                     slope_absolute_error=slope_absolute,
                     slope_relative_l2=slope_relative,
                     reference_raw_slope_l2=float(np.linalg.norm(rs[:, k])),
                     candidate_raw_slope_l2=float(np.linalg.norm(cs[:, k])),
                     reference_flux=rf, candidate_flux=cf,
                     flux_relative_error=flux_error,
                     photon_conservation_relative_error=conservation)
        probes.append(entry)
        checks.extend((gate(name+':rate', rate_error, rate_limit),
                       gate(name+':flux', flux_error, thresholds['flux_relative_error']),
                       gate(name+':slope_absolute', slope_absolute, thresholds['slope_absolute_error'])))
        if name not in ('flat', 'piston_100nm'):
            checks.append(gate(name+':slope_relative', slope_relative,
                               thresholds['excited_slopes_relative_l2']))
        for owner, value in conservation.items():
            checks.append(gate(name+':'+owner+':conservation', value,
                               thresholds['photon_conservation_relative_error']))
    piston = probe_names.index('piston_100nm')
    piston_checks = {}
    for owner, rates in (('reference', rr), ('candidate', cr)):
        value = relative_l2(rates[:, :, piston]-rates[:, :, flat], rates[:, :, flat])
        piston_checks[owner] = value
        checks.append(gate(owner+':piston_invariance', value, thresholds['piston_rate_relative_l2']))
    derived = dict(reference_slopes=rs, candidate_slopes=cs)
    calibration = None
    present = ['calibration_positive' in arrays for arrays in (reference, candidate)]
    if any(present) and not all(present):
        raise ValueError('Both runs must supply central calibration arrays, or neither')
    if all(present):
        matrices = []
        for owner, arrays in (('reference', reference), ('candidate', candidate)):
            positive = readout_stack(arrays['calibration_positive'], divisor)
            negative = readout_stack(arrays['calibration_negative'], divisor)
            matrix = (positive-negative)/(2*poke)
            matrices.append(matrix)
            derived[owner+'_D'] = matrix
        rd, cd = matrices
        difference = cd-rd
        global_error = relative_l2(difference, rd)
        column_norms = np.linalg.norm(rd, axis=0)
        error_norms = np.linalg.norm(difference, axis=0)
        floor = float(thresholds['column_denominator_floor_fraction']*np.max(column_norms))
        if floor <= 0:
            raise ValueError('Reference calibration has no nonzero columns')
        column_relative = error_norms/np.maximum(column_norms, floor)
        weak = np.flatnonzero(column_norms < floor)
        column_limit = thresholds['calibration_column_relative_l2']
        calibration = dict(relative_l2=global_error,
                           reference_column_norms=column_norms.tolist(),
                           column_error_norms=error_norms.tolist(),
                           column_relative_l2=column_relative.tolist(),
                           maximum_column_relative_l2=float(np.max(column_relative)),
                           worst_column_index_zero_based=int(np.argmax(column_relative)),
                           denominator_floor=floor,
                           denominator_floor_fraction=thresholds['column_denominator_floor_fraction'],
                           weak_column_indices_zero_based=weak.tolist(),
                           weak_column_absolute_error_limit=float(column_limit*floor),
                           column_rule='norm(error[:,j]) <= column_limit * max(norm(reference[:,j]), denominator_floor)',
                           central_poke_m=float(poke), precision='Float64')
        global_limit = (thresholds['cpu_cuda_calibration_relative_l2'] if cpu_cuda
                        else thresholds['calibration_relative_l2'])
        checks.extend((gate('calibration:global', global_error, global_limit),
                       gate('calibration:columns', float(np.max(column_relative)), column_limit)))
    return dict(passed=all(check['passed'] for check in checks),
                full_validation=calibration is not None,
                flat_subtraction=False, divisor=float(divisor),
                analytic_incident_rate=float(incident), probes=probes,
                live_input=live_input,
                piston_invariance=piston_checks, calibration=calibration,
                checks=checks), derived


def self_check(thresholds):
    """Discriminating corruption checks on deterministic nonsymmetric data."""
    rng = np.random.default_rng(823)
    names = ['flat', 'piston_100nm', 'tilt_x_100nm']+['excited_'+str(k) for k in range(14)]
    rate = rng.uniform(1, 2, (64, 64, 17)).astype(np.float32)
    rate /= np.sum(rate, axis=(0, 1), keepdims=True)
    rate[:, :, 1] = rate[:, :, 0]
    positive = rng.uniform(1, 2, (64, 64, 277)).astype(np.float32)
    negative = rng.uniform(1, 2, (64, 64, 277)).astype(np.float32)
    reference = dict(rates=rate, calibration_positive=positive, calibration_negative=negative)
    result, _ = compare_arrays(reference, reference, names, 1, 1, 2e-8, thresholds)
    if not result['passed']:
        raise AssertionError('Identical synthetic evidence failed')
    corruptions = dict(gain={**reference, 'rates': rate*1.1},
                       sign_pupil_exchange={**reference, 'rates': rate[::-1, :, :]},
                       transpose={**reference, 'rates': rate.transpose(1, 0, 2)},
                       calibration_sign={**reference, 'calibration_positive': negative,
                                         'calibration_negative': positive})
    output = {}
    for name, candidate in corruptions.items():
        result, _ = compare_arrays(reference, candidate, names, 1, 1, 2e-8, thresholds)
        failures = [c['name'] for c in result['checks'] if not c['passed']]
        if not failures:
            raise AssertionError('Undetected synthetic corruption: '+name)
        output[name] = dict(detected=True, failed_checks=failures)
    stale = {**reference, 'rates': np.repeat(rate[:, :, :1], 17, axis=2)}
    result, _ = compare_arrays(stale, stale, names, 1, 1, 2e-8, thresholds)
    failures = [c['name'] for c in result['checks'] if not c['passed']]
    if result['passed'] or not all(owner+':live_input' in failures for owner in ('reference', 'candidate')):
        raise AssertionError('Identically stale reference/candidate evidence escaped detection')
    output['both_stale'] = dict(detected=True, failed_checks=failures)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--candidate', type=Path)
    parser.add_argument('--fixture', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--cpu-cuda', action='store_true')
    parser.add_argument('--self-check', action='store_true')
    args = parser.parse_args()
    protocol_path = HERE/'protocol.toml'
    protocol = tomllib.loads(protocol_path.read_text())
    if args.self_check:
        print(json.dumps(self_check(protocol['validation']), indent=2))
        return 0
    if any(value is None for value in (args.reference, args.candidate, args.fixture, args.output)):
        parser.error('--reference, --candidate, --fixture and --output are required')
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('Output must be empty; preserve earlier evaluation evidence')
    fixture, arrays = load_fixture(args.fixture)
    fh, ph = digest(args.fixture/'fixture.toml'), digest(protocol_path)
    reference_record, reference, reference_provenance = load_run(args.reference, fh, ph, fixture['aos_revision'])
    candidate_record, candidate, candidate_provenance = load_run(args.candidate, fh, ph, fixture['aos_revision'])
    if args.cpu_cuda:
        if reference_record['engine'] != candidate_record['engine'] or {reference_record['backend'], candidate_record['backend']} != {'cpu', 'cuda'}:
            raise ValueError('--cpu-cuda requires CPU and CUDA records from the same engine')
    if reference_record['prescription'] != candidate_record['prescription']:
        raise ValueError('Different prescriptions must not be evaluated as equivalent')
    contract = fixture['contract']
    divisor = incident_divisor(contract, arrays['amplitude'])
    from common import support_indices
    incident = divisor*len(support_indices()[0])
    result, derived = compare_arrays(reference, candidate, fixture['probe_names'],
                                    divisor, incident, contract['calibration_poke_m'],
                                    protocol['validation'], args.cpu_cuda)
    result.update(reference=reference_provenance, candidate=candidate_provenance,
                  fixture_sha256=fh, protocol_sha256=ph,
                  evaluator_sha256=digest(__file__), cpu_cuda=args.cpu_cuda,
                  readout_sha256=digest(HERE/'common.py'),
                  thresholds=protocol['validation'],
                  scope='raw common I4Q slopes and central calibration; no fitted correction',
                  timing_eligible=result['passed'] and result['full_validation'] and
                  reference_record['prescription'] == 'common')
    args.output.mkdir(parents=True, exist_ok=True)
    result['derived_arrays'] = {name: save_f64(args.output, name, value)
                               for name, value in derived.items()}
    write_json(args.output/'summary.json', result)
    print(json.dumps(dict(passed=result['passed'], full_validation=result['full_validation'],
                          timing_eligible=result['timing_eligible'],
                          failed_checks=[c for c in result['checks'] if not c['passed']],
                          summary=str(args.output/'summary.json')), indent=2))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
