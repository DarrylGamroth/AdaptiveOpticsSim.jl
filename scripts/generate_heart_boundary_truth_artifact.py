#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parents[1]
REVOLT = ROOT.parent / 'REVOLT'
OUTDIR = ROOT / 'benchmarks' / 'results' / 'truth'
OUTFILE = OUTDIR / '2026-04-09-heart-boundary-truth.toml'
MANIFEST = OUTDIR / 'manifest.toml'

COMMON_TOML = REVOLT / 'Julia' / 'contracts' / 'revolt' / 'common.toml'
SHWFS_TOML = REVOLT / 'Julia' / 'contracts' / 'revolt' / 'shwfs.toml'
CAMERAS_TOML = REVOLT / 'Julia' / 'contracts' / 'revolt' / 'cameras.toml'
DM_MAP = REVOLT / 'Julia' / 'assets' / 'revolt_like' / 'revolt_like_dmActuatorMap_277.csv'
SPECULA_YML = REVOLT / 'Python' / 'specula' / 'params_revolt_modal.yml'
OOPAO_SH = REVOLT / 'Python' / 'OOPAO_model' / 'REVOLTtwin_idealSHWFS.py'
OOPAO_PARAM = REVOLT / 'Python' / 'OOPAO_model' / 'parameter_files' / 'parameterFile_REVOLTtwin_ideal.py'
DAO_PUPIL = REVOLT / 'Data' / 'DAO1.2m-pupil-median-stack.tiff'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def parse_toml(path: Path):
    with path.open('rb') as f:
        return tomllib.load(f)


def parse_specula_value(text: str, key: str) -> str:
    m = re.search(rf'^\s*{re.escape(key)}\s*:\s*(.+?)\s*(?:#.*)?$', text, re.MULTILINE)
    if not m:
        raise RuntimeError(f'missing key {key} in {SPECULA_YML}')
    return m.group(1).strip()


def parse_python_assignment(text: str, name: str) -> str:
    m = re.search(rf'{re.escape(name)}\s*=\s*([^#\n]+)', text)
    if not m:
        raise RuntimeError(f'missing assignment {name}')
    return m.group(1).strip()


def parse_param_value(text: str, key: str) -> str:
    m = re.search(rf"param\['{re.escape(key)}'\s*\]\s*=\s*([^#\n]+)", text)
    if not m:
        raise RuntimeError(f'missing parameter {key}')
    return m.group(1).strip()


def parse_dm_map(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append([int(x) for x in line.replace(',', ' ').split()])
    header = rows[0]
    grid = rows[1:]
    vals = [x for row in grid for x in row]
    return {
        'header': header,
        'rows': len(grid),
        'cols': len(grid[0]),
        'nonzero_count': sum(1 for x in vals if x > 0),
        'max_index': max(vals),
    }


def file_description(path: Path) -> str:
    out = subprocess.check_output(['file', str(path)], text=True).strip()
    return out


def toml_scalar(v):
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)
    s = str(v).replace('\\', '\\\\').replace('"', '\\"')
    return f'"{s}"'


def toml_array(seq):
    return '[' + ', '.join(toml_scalar(x) for x in seq) + ']'


def dumps_toml(data: dict) -> str:
    lines: list[str] = []
    for key in ['artifact_id', 'generated_on']:
        lines.append(f'{key} = {toml_scalar(data[key])}')
    lines.append('')
    lines.append('[scope]')
    for k, v in data['scope'].items():
        if isinstance(v, list):
            lines.append(f'{k} = {toml_array(v)}')
        else:
            lines.append(f'{k} = {toml_scalar(v)}')
    lines.append('')
    lines.append('[boundary_truth]')
    for k, v in data['boundary_truth'].items():
        if isinstance(v, list):
            lines.append(f'{k} = {toml_array(v)}')
        else:
            lines.append(f'{k} = {toml_scalar(v)}')
    lines.append('')
    lines.append('[known_differences]')
    for k, v in data['known_differences'].items():
        lines.append(f'{k} = {toml_scalar(v)}')
    lines.append('')
    lines.append('[verdict]')
    for k, v in data['verdict'].items():
        lines.append(f'{k} = {toml_scalar(v)}')
    for name, info in data['provenance'].items():
        lines.append('')
        lines.append(f'[provenance.{name}]')
        for k, v in info.items():
            lines.append(f'{k} = {toml_scalar(v)}')
    return '\n'.join(lines) + '\n'


def main() -> None:
    common = parse_toml(COMMON_TOML)
    shwfs = parse_toml(SHWFS_TOML)
    cameras = parse_toml(CAMERAS_TOML)
    specula_text = SPECULA_YML.read_text()
    oopao_sh_text = OOPAO_SH.read_text()
    oopao_param_text = OOPAO_PARAM.read_text()
    dm_info = parse_dm_map(DM_MAP)

    n_subap_julia = int(shwfs['wfs']['n_subap'])
    n_pix_subap_julia = int(shwfs['wfs']['n_pix_subap'])
    detector_shape = [n_subap_julia * n_pix_subap_julia, n_subap_julia * n_pix_subap_julia]
    light_threshold = float(parse_param_value(oopao_param_text, 'lightThreshold'))

    data = {
        'artifact_id': 'HEART-TRUTH-2026-04-09',
        'generated_on': '2026-04-09',
        'scope': {
            'artifact_kind': 'scientist_owned_heart_boundary_truth',
            'interpretation': 'Scientist-owned REVOLT/DAO boundary and geometry sources used as external truth for the maintained HEART RTC HIL contract',
            'truth_level': 'boundary_and_geometry',
            'raw_on_sky_telemetry_included': False,
            'provenance_root': str(REVOLT),
        },
        'boundary_truth': {
            'external_command_length': int(dm_info['max_index']),
            'dm_grid_shape': [int(dm_info['rows']), int(dm_info['cols'])],
            'dm_nonzero_actuators': int(dm_info['nonzero_count']),
            'telescope_diameter_m': float(common['telescope']['diameter']),
            'sh_n_subap': n_subap_julia,
            'sh_pixels_per_subap': n_pix_subap_julia,
            'detector_output_shape': detector_shape,
            'detector_response_model': 'cblue1_cmos',
            'julia_threshold': float(shwfs['wfs']['threshold']),
            'oopao_light_threshold': light_threshold,
            'specula_subap_wanted_fov_arcsec': float(parse_specula_value(specula_text, 'subap_wanted_fov')),
            'specula_sensor_pxscale_arcsec': float(parse_specula_value(specula_text, 'sensor_pxscale')),
            'specula_detector_size': [352, 352],
        },
        'known_differences': {
            'central_obstruction_note': 'Julia common contract uses 0.0, while SPECULA scientist-owned config records an approximate 0.25 m obstruction for DAO; this remains a known model-normalization difference rather than a boundary mismatch.',
            'wfs_wavelength_note': 'SPECULA SH config is explicit at 800 nm; OOPAO SH path uses R-band source and Julia uses REVOLT visible-star plus cblue1 camera settings. This artifact treats command/frame boundary and geometry as truth, not full optical identity.',
            'oopao_sh_construction_note': 'Scientist-owned OOPAO SH script sets n_subaperture=16 and shannon_sampling=true for the ideal SHWFS path.',
        },
        'verdict': {
            'boundary_truth_available': True,
            'suitable_for_supported_boundary_claims': True,
            'suitable_for_full_optical_truth_claims': False,
        },
        'provenance': {
            'julia_common_toml': {
                'path': str(COMMON_TOML),
                'sha256': sha256(COMMON_TOML),
            },
            'julia_shwfs_toml': {
                'path': str(SHWFS_TOML),
                'sha256': sha256(SHWFS_TOML),
            },
            'julia_cameras_toml': {
                'path': str(CAMERAS_TOML),
                'sha256': sha256(CAMERAS_TOML),
            },
            'dm_map': {
                'path': str(DM_MAP),
                'sha256': sha256(DM_MAP),
            },
            'specula_params': {
                'path': str(SPECULA_YML),
                'sha256': sha256(SPECULA_YML),
            },
            'oopao_sh_script': {
                'path': str(OOPAO_SH),
                'sha256': sha256(OOPAO_SH),
            },
            'oopao_parameter_file': {
                'path': str(OOPAO_PARAM),
                'sha256': sha256(OOPAO_PARAM),
            },
            'dao_pupil_image': {
                'path': str(DAO_PUPIL),
                'sha256': sha256(DAO_PUPIL),
                'file_description': file_description(DAO_PUPIL),
            },
        },
    }

    OUTDIR.mkdir(parents=True, exist_ok=True)
    OUTFILE.write_text(dumps_toml(data))

    manifest = {'artifacts': []}
    if MANIFEST.exists():
        import tomllib as _tomllib
        manifest = _tomllib.loads(MANIFEST.read_text())
    artifacts = [a for a in manifest.get('artifacts', []) if a.get('id') != 'HEART-TRUTH-2026-04-09']
    artifacts.append({
        'purpose': 'scientist-owned HEART boundary truth artifact from REVOLT/DAO assets',
        'id': 'HEART-TRUTH-2026-04-09',
        'path': OUTFILE.name,
    })
    MANIFEST.write_text('[[artifacts]]\n' + '\n\n[[artifacts]]\n'.join(
        '\n'.join(f'{k} = {json.dumps(v)}' for k, v in art.items()) for art in artifacts
    ) + '\n')
    print(OUTFILE)


if __name__ == '__main__':
    main()
