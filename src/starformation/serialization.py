"""Versioned, exclusive scientific datasets with a completion marker."""

import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from .cloud import json_ready
from .engine import SimulationResult
from .visualization import validate_dataset

SCHEMA_VERSION = 1
_FILES = ('configuration.json', 'initial_state.npz', 'trajectory.npz',
          'diagnostics.json', 'validation.json', 'metadata.json')


def _write_json(path, payload):
    with path.open('x') as stream:
        json.dump(json_ready(payload), stream, indent=2, allow_nan=False)
        stream.write('\n')


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _revision():
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=root,
                            capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def save_result(result: SimulationResult, output):
    """Never overwrite. A directory lacking COMPLETE.json is not a dataset."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    config = result.config
    _write_json(output/'configuration.json', {
        'schema_version': 1, 'units': 'SI', 'parameters': config.as_si_dict(),
        'input_units': {'cloud_mass_msun': 'M_sun', 'radius_pc': 'pc',
                        'dt_over_initial_tff': 'dimensionless',
                        'softening_over_radius': 'dimensionless'}})
    np.savez_compressed(output/'initial_state.npz', pos=result.initial.pos,
                        vel=result.initial.vel, mass=result.initial.mass)
    np.savez_compressed(output/'trajectory.npz', **result.trajectory)
    report = dict(result.report)
    report['trajectory_sha256'] = _sha(output/'trajectory.npz')
    _write_json(output/'diagnostics.json', report)
    accuracy = result.report['accuracy']
    _write_json(output/'validation.json', {
        'schema_version': 1, 'criterion': 'peak sampled |delta E| <= energy_budget',
        'status': accuracy['energy_status'],
        'peak_energy_error': accuracy['peak_energy_error'],
        'energy_budget': accuracy['energy_budget'],
        'sample_start_s': accuracy['sample_start_s'],
        'sample_end_s': accuracy['sample_end_s'],
        'warnings': (['Energy criterion FAILED on sampled interval.']
                     if accuracy['energy_status'] == 'FAIL' else []) +
                    ['Legacy dt=0.012 initial t_ff failed the Phase 1 1% criterion over 1.2 t_ff.',
                     'No full-duration or hydrodynamic validation is implied.']})
    _write_json(output/'metadata.json', {
        'schema_version': 1, 'code_revision': _revision(),
        'requested_output_path': str(output.resolve()),
        'python': sys.version.split()[0],
        'dependencies': {p: importlib.metadata.version(p) for p in
                         ('numpy', 'astropy', 'matplotlib')},
        'seed': result.initial.seed, 'initialization_method': result.initial.method,
        'softening': {'law': 'Plummer', 'length_m': config.softening_m},
        'timestep_history': {'kind': 'fixed', 'dt_s': config.dt_s,
                             'completed_steps': config.steps},
        'validation_status': accuracy['energy_status']})
    _write_json(output/'COMPLETE.json', {
        'schema_version': 1, 'status': 'complete',
        'artifact_schemas': {'configuration': 1, 'initial_state': 1,
                             'particle_snapshots': 1, 'diagnostics': 1,
                             'validation': 1, 'metadata': 1},
        'files_sha256': {name: _sha(output/name) for name in _FILES}})
    return output


def load_result(directory):
    """Load trusted numeric arrays only; reject incomplete or altered datasets."""
    directory = Path(directory)
    marker_path = directory/'COMPLETE.json'
    if not marker_path.is_file():
        raise ValueError('Incomplete dataset: COMPLETE.json is absent')
    marker = json.loads(marker_path.read_text())
    if marker.get('schema_version') != SCHEMA_VERSION or marker.get('status') != 'complete':
        raise ValueError('Unsupported or incomplete dataset')
    if marker.get('artifact_schemas') != {'configuration': 1, 'initial_state': 1,
                                          'particle_snapshots': 1, 'diagnostics': 1,
                                          'validation': 1, 'metadata': 1}:
        raise ValueError('Unsupported artifact schema version')
    if set(marker.get('files_sha256', {})) != set(_FILES):
        raise ValueError('Dataset file manifest mismatch')
    for name in _FILES:
        if _sha(directory/name) != marker['files_sha256'][name]:
            raise ValueError(f'Dataset checksum mismatch: {name}')
    config = json.loads((directory/'configuration.json').read_text())
    report = json.loads((directory/'diagnostics.json').read_text())
    validation = json.loads((directory/'validation.json').read_text())
    metadata = json.loads((directory/'metadata.json').read_text())
    for name, value in [('configuration', config), ('diagnostics', report),
                        ('validation', validation), ('metadata', metadata)]:
        if value.get('schema_version') != SCHEMA_VERSION:
            raise ValueError(f'Unsupported {name} schema version')
    if config.get('units') != 'SI' or report.get('units') != 'SI':
        raise ValueError('Only SI datasets are supported')
    with np.load(directory/'initial_state.npz', allow_pickle=False) as stored:
        initial = {k: stored[k] for k in stored.files}
    with np.load(directory/'trajectory.npz', allow_pickle=False) as stored:
        data = {k: stored[k] for k in stored.files}
    if set(initial) != {'pos', 'vel', 'mass'} or set(data) != {
            'pos', 'vel', 'mass', 'time_s', 'density_kg_m3',
            'neighbour_radius_m', 'radial_velocity_m_s'}:
        raise ValueError('Unexpected array fields')
    for key in initial:
        if not np.array_equal(initial[key], data[key][0] if key != 'mass' else data['mass']):
            raise ValueError(f'Initial state mismatch: {key}')
    validate_dataset(data, report)
    if validation['status'] != report['accuracy']['energy_status'] or metadata['validation_status'] != validation['status']:
        raise ValueError('Validation status mismatch')
    if config['parameters']['dt_s'] != report['configuration']['dt_s']:
        raise ValueError('Configuration timestep mismatch')
    return {'configuration': config, 'initial_state': initial, 'trajectory': data,
            'diagnostics': report, 'validation': validation, 'metadata': metadata}
