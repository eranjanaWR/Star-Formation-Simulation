"""Render an already completed saved dataset; never launches dynamics."""

import argparse
from pathlib import Path

from .serialization import load_result
from .visualization import render


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--view', choices=('com', 'fixed'))
    parser.add_argument('--animation', action='store_true')
    args = parser.parse_args(argv)
    saved = load_result(args.data)
    configured = saved['configuration']['parameters']['render_path']
    if configured is not None and Path(configured) != args.output:
        parser.error('--output differs from configured output.render_path')
    settings = saved['configuration']['parameters']
    render(args.data, args.output, view=args.view or settings['render_view'],
           animation=args.animation or settings['render_animation'],
           fps=settings['render_fps'], dpi=settings['render_dpi'])
    print(f'{args.output}: rendered {saved["validation"]["status"]} dataset')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
