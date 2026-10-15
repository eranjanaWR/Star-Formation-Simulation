"""Safe source-tree bridge to the packaged code.

The preserved ``starformation.py`` has import-time execution. Python gives this
package directory precedence over that file, so source-tree imports stay safe.
Installed distributions use ``src/starformation`` directly.
"""

from pathlib import Path

__path__ = [str(Path(__file__).resolve().parents[1] / 'src' / 'starformation')]
__version__ = '0.3.0'
