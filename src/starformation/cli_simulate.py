"""Explicit gravity-only simulation CLI; never imports Matplotlib."""

import argparse
from pathlib import Path

from .config import load_config
from .engine import run
from .serialization import save_result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    config = load_config(args.config)
    if config.simulation_path is not None and Path(config.simulation_path) != args.output:
        parser.error('--output differs from configured output.simulation_path')
    if args.output.exists():
        parser.error(f'Output already exists: {args.output}')
    result = run(config)
    save_result(result, args.output)
    accuracy = result.report['accuracy']
    print(f'{args.output}: {accuracy["energy_status"]}; peak sampled energy error '
          f'{accuracy["peak_energy_error"]:.6%}')
    return 1 if accuracy['energy_status'] == 'FAIL' else 0


if __name__ == '__main__':
    raise SystemExit(main())
