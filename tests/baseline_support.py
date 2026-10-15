"""Load only the frozen baseline's initialization and three physics functions.

The legacy module cannot be imported safely: import runs the full simulation.
This test-only adapter executes an AST prefix from an exact, hash-pinned source.
It deliberately excludes plotting imports, output, staging, and rendering.
"""

import ast
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "baseline" / "manifest.json"


def load_baseline():
    manifest = json.loads(MANIFEST.read_text())
    path = ROOT / "starformation.py"
    source = path.read_bytes()
    expected = manifest["preserved_files"]["starformation.py"]["sha256"]
    if hashlib.sha256(source).hexdigest() != expected:
        raise RuntimeError("Baseline source changed; review the adapter before executing it")
    tree = ast.parse(source, filename=str(path))
    selected = []
    found_boundary = False
    for node in tree.body:
        if isinstance(node, ast.Import):
            if all(alias.name == "numpy" for alias in node.names):
                selected.append(node)
        elif isinstance(node, ast.ImportFrom):
            if node.module == "astropy":
                selected.append(node)
        elif isinstance(node, (ast.Assign, ast.FunctionDef)):
            selected.append(node)
        elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            # Preserve the actual np.random.seed(42) call, not just assignments.
            if ast.unparse(node.value.func) == "np.random.seed":
                selected.append(node)
        if isinstance(node, ast.FunctionDef) and node.name == "compute_jeans_length":
            found_boundary = True
            break
    if not found_boundary:
        raise RuntimeError("Baseline extraction boundary missing")
    state = np.random.get_state()
    namespace = {"__name__": "baseline_under_test"}
    try:
        exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), "exec"), namespace)
    finally:
        np.random.set_state(state)
    return namespace


def summary(namespace):
    """Small numerical fingerprint, not a conservation or convergence benchmark."""
    n = namespace
    p, v, m = n["pos"], n["vel"], n["mass"]
    a = n["compute_acceleration"](p, m, n["SOFTENING"])
    p1, v1, a1 = n["leapfrog_step"](p, v, a, n["dt"])
    au = n["u"].au.to(n["u"].m)
    return {
        "initial_position_first3_au": (p[:3] / au).tolist(),
        "initial_velocity_first3_m_s": v[:3].tolist(),
        "initial_com_au": (np.sum(m[:, None] * p, axis=0) / np.sum(m) / au).tolist(),
        "initial_com_velocity_m_s": (np.sum(m[:, None] * v, axis=0) / np.sum(m)).tolist(),
        "initial_position_rms_au": float(np.sqrt(np.mean(p * p)) / au),
        "initial_speed_rms_m_s": float(np.sqrt(np.mean(np.sum(v * v, axis=1)))),
        "initial_legacy_jeans_au": float(n["compute_jeans_length"](p, m) / au),
        "first_step_position_first3_au": (p1[:3] / au).tolist(),
        "first_step_velocity_first3_m_s": v1[:3].tolist(),
        "first_step_acceleration_first3_m_s2": a1[:3].tolist(),
    }
