# Phase 3 — modular scientific research architecture

Phase 3 restructures the **same** Plummer-softened, collisionless N-body model.
It adds no gas dynamics, pressure, thermal evolution, sinks, or protostar claims.
It does not improve the accuracy of the force law or timestep. The Phase 1
legacy-step failure remains a failure.

## Preflight and migration

Work began on branch `dev` at committed Phase 2 revision `92986ac` with a clean
working tree. The Phase 0–2 reports, corrected seed-42 fixture, source/media
manifest, archived Phase 1/2 results, dependencies, and tests were inspected.
All 53 inherited tests passed before refactoring. The preserved script and both
animation hashes are unchanged, as enforced by the Phase 0 tests.

Dependency direction (arrows mean imports or calls):

```text
config -> units, initial_conditions -> config
gravity/integration/energy -> _nbody_impl
cloud -> gravity + energy
engine -> config + initial_conditions + gravity + integration + cloud
serialization -> engine type + cloud JSON converter + visualization validator
cli_simulate -> config + engine + serialization
cli_render -> serialization + visualization
visualization -> NumPy; Matplotlib only inside plotting functions
root compatibility modules -> packaged modules
```

The `src/starformation` package contains independent scientific modules. The
root `nbody.py`, `scientific_diagnostics.py`, and `scientific_visualization.py`
are compatibility wrappers. The historical `starformation.py` remains untouched.
Because that script has import-time execution and shares the package's name, the
small root `starformation/__init__.py` source-tree bridge takes import precedence;
an installed distribution uses `src/starformation` directly. Explicitly running
`python starformation.py` still invokes the hazardous historical full job and is
not a supported research CLI.

| Previous component | New location | Compatibility decision |
| --- | --- | --- |
| `nbody.py` force, KDK, energy | `_nbody_impl.py`, `gravity.py`, `integration.py`, `energy.py` | Exact arithmetic body copied; root signatures retained |
| `scientific_diagnostics.py` | `cloud.py` | Only import changed to package-relative; root functions re-exported |
| `scientific_visualization.py` | `visualization.py` | Plot code copied; root CLI/function names retained |
| Phase 0 scripted initialization | `initial_conditions.py` | Draw order and formulas reproduced with isolated legacy MT19937 |
| Phase 2 preview orchestration | `engine.py`, `serialization.py`, `cli_simulate.py` | General controlled fixed-step runs; archived preview remains unchanged |
| Phase 2 data renderer | `cli_render.py` | Reads completed data only |

Compatibility risks were RNG sequence, floating-point expression order, source
import collision, plotting side effects, changed public imports, and accidentally
relabelled numerical failures. They have direct regression tests. Historical
scripts and data are not migrated in place; the Phase 2 dataset format remains
readable through the old visualization API.

## Public interfaces and units

The force, KDK, potential, and conservation functions retain their Phase 1
signatures: `compute_acceleration(pos,mass,softening,*,G)`,
`leapfrog_step(pos,vel,acc,dt,mass,softening,*,G)`,
`potential_energy(pos,mass,softening,*,G)`, and
`diagnostics(pos,vel,mass,softening,*,G)`. Their numeric inputs are SI and their
validated arithmetic is unchanged. `cloud.snapshot(...)` and
`cloud.accuracy_history(...)` retain Phase 2 definitions and 1% default energy
budget. `make_initial_state(config, *, pos=None, vel=None, mass=None)` returns a
state and method/seed. All three explicit arrays must be provided together;
particle masses must sum to the configured mass. COM position/velocity subtraction
is opt-in, never a default correction. `engine.run(config, initial=None)` returns
every state from t=0 through the requested final step without writing or plotting.
`save_result(result, output)` writes an exclusive new directory;
`load_result(directory)` validates and loads it.

The JSON boundary uses solar masses and parsecs; `units.py` uses the installed
Astropy conversions to SI. Inputs name each unit rather than accepting implicit
unit-bearing solver arrays. Internally positions/softening/radius are m, masses
kg, velocities m/s, times s, acceleration m/s², G m³/(kg s²), energy J, and
momenta SI. The legacy preset retains N=300, seed=42, 1 solar mass, 0.1 pc,
softening=0.05R, dt=0.012 initial t_ff, 1080 steps, Gaussian sigma=0.05
sqrt(GM/R), and rotation=0.03 sqrt(GM/R³). It is a reference preset, **not**
an automatically launched run or a validated timestep. The explicit short-run
config instead uses dt=0.006 t_ff and eight steps; its short duration cannot
claim the 1.2 t_ff Phase 1 criterion.
The new engine saves both t=0 and the state after its last requested step:
`steps=1080` therefore means 1081 saved states. The historical movie stored
1080 pre-advance frames and computed an additional, unstored last state; its
frame convention is preserved only in the historical script and archives.

`config.py` strictly validates the version-1 JSON object and nested `output` and
`render` objects. Required input keys are `schema_version`, `cloud_mass_msun`,
`radius_pc`, `particle_count`, `seed`, `gravitational_constant_si`,
`softening_over_radius`, `dt_over_initial_tff`, `steps`,
`velocity_sigma_over_virial`, `rotation_omega_over_dynamical`, `recenter_com`,
`recenter_velocity`, `sound_speed_m_s`, `neighbours`, `energy_budget`, `output`,
and `render`. `output` contains `simulation_path` and `render_path` (null or an
explicit path); command-line `--output` is always required and must agree if a
path is configured. `render` has `view` (`com`/`fixed`), `animation`, `fps`, and
`dpi`. The new render CLI uses these defaults, with `--view` and `--animation`
overrides; fps and dpi control optional preview and figure output. The seed is a uint32
MT19937 seed; steps, counts, k, physical values, and booleans are range-checked.
Unknown fields or schema versions are rejected.

## Versioned dataset and provenance

Every new experiment directory is exclusive. `COMPLETE.json` is written **last**
and declares dataset/artifact schema version 1 plus SHA-256 for every required
file. A loader rejects a missing marker, unknown version, hash mismatch,
malformed numeric state, inconsistent initial state, inconsistent energy history,
or changed PASS/FAIL label. No pickle is loaded (`allow_pickle=False`). Arrays
remain binary NumPy float64 through compressed NPZ, without JSON rounding.

| File | Contents |
| --- | --- |
| `configuration.json` | version, SI parameter record, input-unit map |
| `initial_state.npz` | initial `(N,3)` position/velocity and `(N,)` mass |
| `trajectory.npz` | `(steps+1,N,3)` states; `(steps+1,)` time; particle density/support/radial flow |
| `diagnostics.json` | Phase 2-compatible sampled summaries, definitions, energy/residual histories, trajectory checksum |
| `validation.json` | 1% or configured energy criterion, sampled interval, status and warnings |
| `metadata.json` | Git revision, Python/direct dependency versions, seed/method, Plummer softening, fixed timestep history, requested output path, status |
| `COMPLETE.json` | completion status, artifact schema map, file checksums |

The `diagnostics.json` energy error is `(E-E0)/|E0|`; the running peak is the
maximum absolute error over saved samples. Momentum, angular momentum, and COM
residuals retain Phase 1/2 normalizations. A `FAIL` status is saved rather than
adjusting energy or momentum. The old Phase 2 archives are schema version 1 in
their own format and remain read-only; the new loader intentionally requires the
new completion manifest rather than silently guessing a format. Future changes
need a new schema version and an explicit migration reader, not rewriting old
archives. A process interruption can leave a partial directory, but no marker;
such a directory is not loadable or silently overwritten.

## Reproduction and execution

Tested on Python 3.13.5 with NumPy 2.4.4, Astropy 7.2.0, Matplotlib 3.10.8.
`pyproject.toml` supports Python `>=3.13,<3.14`, NumPy `>=2.4,<3`, Astropy
`>=7.2,<8`; rendering has optional Matplotlib `>=3.10,<4`. This is a tested
environment and supported-range declaration, not a clean-environment or
cross-platform certification. Installation, if desired, is
`.venv/bin/python -m pip install -e '.[render]'`; no installation was needed for
the verified source-tree commands:

```sh
.venv/bin/python -B -m unittest discover -s tests -v
.venv/bin/python -B -m starformation.cli_simulate \
  --config configs/phase3_short.json --output /private/tmp/my-new-run
.venv/bin/python -B -m starformation.cli_render /private/tmp/my-new-run \
  --output /private/tmp/my-new-figures --view fixed
```

Both output paths must be new; rendering cannot initiate a simulation and
simulation does not import Matplotlib or create media. `--animation` on the
render CLI opts into the existing short preview GIF, in a new output directory.
The older `scientific_visualization.py` command still renders Phase 2 archives.
For the full existing numerical references, inspect `docs/PHASE1_NBODY_VALIDATION.md`,
`docs/PHASE2_SCIENTIFIC_DIAGNOSTICS.md`, and `validation/phase1`,
`validation/phase2`. No full 1080-step run was performed in Phase 3.

## Verification results

Acceptance rules were set before checks: seed-42 initial arrays and the archived
100-step trajectory must match bitwise; diagnostic energy and peak error must
match their recorded values within floating-point representation; all 53
inherited tests and new tests must pass; source/media hashes must match the
Phase 0 manifest; simulation import must not import Matplotlib. The previous
1% peak sampled-energy budget remains the numerical criterion, not a physical
model-validation certificate.

- Seed-42 initial position, velocity, and mass arrays: bitwise equal to
  `validation/phase1/cloud_initial.npz`; global NumPy RNG remains unchanged.
- Legacy dt=0.012 t_ff, 100-step trajectory: position, velocity, mass, and time
  arrays bitwise equal to `validation/phase2/dt0.012/trajectory.npz`; sampled
  kinetic/potential/total energy, r50, and inflow fraction match exactly.
- Its peak sampled energy error is **0.019005008127973757 (1.9005008%)**, still
  `FAIL` at the 1% budget; serialized status remains `FAIL`.
- An explicit eight-step seed-42 dt=0.006 t_ff smoke run completed, recording
  nine states and **0.000005%** peak sampled energy error (rounded CLI value).
  Static fixed-origin figures were generated separately from its saved data in
  a temporary directory. This very short pass says nothing about later collapse.
- All 53 inherited tests and 10 new architecture tests pass (63 total); the
  inherited suite was also run before any modifications. Import isolation,
  deterministic replay, strict configuration, exclusive paths, version rejection,
  serialization round-trip, and archived numerical equivalence are covered.
- `starformation.py`, MP4, and GIF SHA-256 values remain respectively
  `a30fd78d36c72141d1716a08e6c06dd0f2b80d658a0b785ebe39d1012ea6f438`,
  `78baa44c73e36f64b4bc441cfea38977e870fdcd69858d2cbf16dcb8464888c2`,
  and `27e3d5dc267f51d3e02dcbe87bf8f900c3b4e9bbddb88bc7625fc7eb56b5d1b5`.

No Phase 3 test failed after implementation. The underlying **scientific** failure
is the unchanged legacy timestep's 1.9% energy excursion over the tested
1.2 t_ff. The short dt=0.006 pass does not validate a full trajectory.

## Limitations and Phase 4 boundary

The model is collisionless particles with softened gravity; neighbour density
and Jeans values remain diagnostic proxies, not evolved gas. The direct solver
is O(N²); fixed dt may underresolve close encounters; neither full legacy
duration nor cross-platform bitwise reproducibility has been demonstrated.
Astropy constant definitions and dependency versions should be recorded for each
run. The root Phase 2 CLI retains its original rendering defaults. NPZ plus JSON validates stored
consistency and checksums but cannot prove external data's scientific truth.

Before any separately approved SPH work, define an explicit gas state and unit
contract, force composition separate from this reference solver, independent
conservation tests, resolution/convergence criteria, and a new versioned data
schema. Do not reinterpret these N-body results as evidence of star formation.
