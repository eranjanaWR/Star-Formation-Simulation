# Phase 0 — preserved and characterized baseline

Recorded 2026-10-15. Authorization covers Phase 0 only. This record distinguishes
source preservation, measured regression behavior, and independent physical
references. It does not certify the gravitational solver or star-formation physics.

## Implementation summary

The original `starformation.py`, MP4, and GIF are unchanged. A test-only adapter
loads initialization and three physics functions without importing the executable
script. Standard-library `unittest` provides lightweight tests; no packages were
installed. A manifest and numerical fixture record the reference state. The README
has a prominent status notice; its historical body is retained, including claims
that this record explicitly supersedes. No physics or animation behavior changed.

## Preservation and version-control strategy

The starting tree was clean at commit
`774d69713da47c5578badb30a5981d80319bd67f`. That existing immutable commit contains all
four original tracked files and is the historical reference. No reset, checkout,
history rewrite, automatic commit, branch change, or duplicate media copy was needed.
Phase 0 additions remain ordinary reviewable working-tree changes.

`baseline/manifest.json` records SHA-256 hashes of the source and both media files,
the original README hash, parameters, tolerances, media metadata, and environment.
Tests enforce exact byte preservation of the source and media. The README changes
intentionally, so its historical hash is provenance rather than a current-file test.

The original README can be inspected without modifying the tree:

```sh
git show 774d69713da47c5578badb30a5981d80319bd67f:README.md
```

Commit the reviewed Phase 0 additions through the project's normal version-control
workflow. In a future approved refactor, keep this baseline available from the
recorded commit and explicitly adapt regression loading; do not silently update
checksums or reference numbers to make tests pass.

## Environment and minimal test configuration

Observed Python: 3.13.5. Direct dependencies: NumPy 2.4.4, Astropy 7.2.0,
Matplotlib 3.10.8. FFmpeg 8.1 is installed with libx264; pip is 26.0.1.
`requirements-baseline.txt` captures installed direct and transitive packages.
It is an environment snapshot, not a wheel/hash lock, platform specification, or
proof of the environment used to produce the historical movies. An independent
clean-environment recreation has not been performed.

From the repository root, use the existing environment:

```sh
.venv/bin/python -B -m unittest discover -s tests -v
```

No installation, display, plotting import, encoder, network access, or output
directory is needed for this command. `-B` prevents bytecode writes. The test suite
reads archived media for hashes and runs only one-step smoke/regression evaluations.
The adapter preserves and restores the caller's NumPy RNG state. It checks the
source hash before extracting an AST prefix ending at `compute_jeans_length`;
staging, plotting, and saving are excluded. It includes the actual seed call.
This is a temporary characterization mechanism, not the proposed modular solver.

The JSON manifest and standard unittest discovery command are sufficient minimal
configuration for Phase 0. A packaging system, test-runner dependency, and final
module architecture would add unnecessary scope here.

## Recorded setup and independently calculated quantities

| Quantity | Baseline value | Meaning |
| --- | --- | --- |
| Total mass | 1 nominal solar mass = 1.988409870698051e30 kg | Equal particle masses |
| Radius | 0.1 pc = 3.0856775814913675e15 m = 20626.480625 AU | Initial bounding sphere |
| Particles | 300 | Collisionless gravity particles |
| Particle mass | 6.628032902326836e27 kg | Total mass / 300 |
| Sound speed | 200 m/s | Annotation only; no pressure acceleration |
| Seed | 42 | Legacy NumPy RandomState/MT19937 |
| Gaussian dispersion | 0.05 sqrt(GM/R) per component | Independent random velocities |
| Angular velocity | 0.03 sqrt(GM/R^3) about z | Added solid-body rotation |
| Softening | 0.05 R = 1031.324031 AU | Fixed Plummer softening |
| Initial mean density | 1.615718379859146e-17 kg/m^3 | Total mass / full sphere volume |
| Initial free-fall time | 523685.938654 yr | Pressureless uniform-sphere reference |
| Timestep | 0.012 initial t_ff = 6284.231264 yr | Fixed, not adaptive |
| Frames | 1080 | State stored before each advance |
| Last stored time | 6780685.533687 yr = 12.948 t_ff | Frame index 1079 |
| Last computed, unstored time | 6786969.764950 yr = 12.96 t_ff | Advance after frame 1079 |
| Uniform-density Jeans length | 72159.551535 AU | Independent diagnostic reference |
| Uniform-density Jeans mass | 5.351999816 solar masses | Mass within radius lambda_J/2 |

The independent test uses explicit SI constants rather than the script's Astropy
expressions: CODATA 2018 G = 6.67430e-11 m^3 kg^-1 s^-2, nominal solar GM =
1.3271244e20 m^3 s^-2, AU = 149597870700 m, and Julian year = 31557600 s.
It derives M_sun = GM_sun/G and pc = AU * 648000/pi. These match the installed
Astropy constants and units; their origins and versioning are documented in
[Astropy's constants reference](https://docs.astropy.org/en/stable/constants/index.html).
Numerical tolerances measure floating-point agreement, not physical uncertainty in G.

For the uniform-sphere analytical reference:

```text
rho_0 = M / [(4*pi/3) R^3]
t_ff = sqrt[3*pi/(32*G*rho_0)] = pi/sqrt(8) * sqrt[R^3/(G*M)]
lambda_J = sqrt[pi*c_s^2/(G*rho_0)]
M_J = (4*pi/3) rho_0 (lambda_J/2)^3
```

G*rho has units s^-2, so t_ff is a time and lambda_J is a length. The Jeans mass
convention is stated explicitly. The free-fall reference assumes an initially
static, pressureless, uniform sphere, unlike the rotating, softened discrete
realization. Jeans theory concerns linear perturbations in an idealized uniform
gas: wavelengths larger than lambda_J are unstable. A finite cloud needs boundary
conditions and a suitable stability analysis; a global length alone cannot classify
all local modes. The [FLASH documentation](https://flash.rochester.edu/site/flashcode/user_support/flash_ug_devel/node139.html)
gives the Jeans length and free-fall formulas and explains resolution constraints.
Here c_s does not enter particle evolution, so this is not a dynamical gas-stability test.

## Seeded measurements and earlier audit correction

The actual seed-42 initial centre of mass is (-187.401418, 985.701306,
-171.758637) AU, with mean velocity (1.252939, 0.124499, 0.591672) m/s.
The legacy Jeans estimator returns 59490.883015 AU initially. That is distinct
from the full-sphere reference above because it assigns all mass inside r_68.
Initial three-dimensional RMS speed is 18.255102 m/s.

`baseline/reference.json` records selected initial positions/velocities, global
summaries, and selected positions/velocities/accelerations after one KDK step.
These are observations from the legacy implementation, not an analytical solution.
The fixture samples the state; it is not a complete saved 300-particle trajectory.
All initial arrays are also checked for exact repeatability in the current environment.

The previous conversational audit's extraction harness omitted the seed call when
selecting assignments and function definitions. Its reported realization-dependent
COM, velocities, virial ratio, pair separations, initial Jeans estimate, and
100-step conservation measurements must not be treated as seed-42 baseline results.
This supersedes the earlier suggestion that the different initial Jeans value showed
an artifact reproducibility mismatch: 59490.883 AU matches the archived MP4's rounded
59490.9 AU annotation. Matching one annotation does not establish full movie provenance.
The independent mass/radius/free-fall calculations were unaffected. The prior
64-particle experiment was separate and is not Phase 0 baseline validation.

## Existing calculation and rendering workflow

The script seeds NumPy, defines SI quantities, samples a uniform sphere using
r = R U^(1/3) and isotropic angles, then adds Gaussian velocities and weak rotation.
All masses are equal. The variables L, M, and T are unused; integration remains in
SI despite the comment about normalized units.

The force routine sums Plummer-softened gravity. The integrator uses half kick,
full drift, recomputed acceleration, and half kick with fixed dt and global mass
and softening. Phase 0 preserves those functions without verifying their accuracy.
Positions, velocities, and the legacy Jeans annotation are stored for every frame.
The rendering phase draws fixed-limit XY and XZ projections, stylizes particles
using radius from the origin, adds labels, and encodes an MP4 at 18 fps and 120 dpi.

There is no main guard. Running **or importing** the script triggers all 1080
advances and rendering. `OUT` is relative to the current working directory, so a
normal run from the repository can overwrite the tracked movie. The module
docstring's filename and GIF output instructions are stale. The script writes
only MP4; no GIF-generation workflow is present.

Historical media inspected with ffprobe:

| Artifact | Codec | Dimensions | Frames | Rate | Duration | Bytes |
| --- | --- | --- | --- | --- | --- | --- |
| MP4 | H.264 | 1440 x 720 | 1080 | 18 fps | 60 s | 4401663 |
| GIF | GIF | 1440 x 720 | 120 | 20 fps | 6 s | 3406207 |

These are preserved reference artifacts, not regenerated validation outputs.
Full movie rerendering and pixel/trajectory equivalence remain untested.
If a full run is approved later, execute an absolute path to the source from a
new, empty temporary directory with a noninteractive Matplotlib backend and a
temporary Matplotlib cache. Confirm the working directory first; outputs will be
written there. Do not run it in the repository or import it as a test helper.

## Validation results and predefined criteria

Criteria were stated before implementation: exact preservation hashes; 1e-12
relative tolerance and zero absolute tolerance for finite nonzero numerical
regression values; 1e-10 relative tolerance for independent SI calculations;
correct shapes and finite values. These tolerances allow rounding differences
while detecting meaningful changes in the baseline. They are not accuracy targets
for a gravitational orbit or collapsing cloud. Exact equality is used only for
same-environment seeded repeats and preserved bytes.

On the recorded environment, all 10 tests passed in approximately 0.26 seconds:

| Test | Observed result | Criterion |
| --- | --- | --- |
| Source and two media hashes/sizes | All match | Exact |
| Full source syntax | Parses and compiles | No execution |
| Recorded scalar parameters | Match source | Exact |
| Independent SI constants and t_ff/dt | Max relative error 2.22e-16 | <= 1e-10 |
| Initial arrays, bounds, mass | Finite; expected shapes; mass error 2.22e-16 | Relative mass error <= 1e-12 |
| Seeded repeat and caller RNG state | Identical | Exact |
| Initial/one-step regression fixture | Max relative error 0 | <= 1e-12 |
| One-step smoke and input preservation | Finite (300,3) results; inputs unchanged | Shapes/finiteness/exact inputs |
| Last stored frame and playback duration | 12.948 t_ff and 60 seconds | <= 1e-10 or exact duration |
| Adapter isolation in fresh process | No Matplotlib or staging objects | Exit status zero |

`git diff --check` also passed. No full simulation, long-duration conservation test,
renderer invocation, dependency installation, or package upgrade was performed.
The pip environment inspection emitted a disabled-cache warning because its default
cache is not writable; package listing succeeded and no environment change was made.

## Scientific limitations and next phase

The model has no pressure, hydrodynamics, thermal evolution, cooling, shocks, sinks,
or hydrostatic stellar structure. It cannot establish molecular-cloud fragmentation,
accretion disks, or protostar birth. Uniform sampling is not a pressure-confined
Bonnor-Ebert equilibrium. Softening is not an SPH smoothing length.

The fixed-step KDK scheme is not exactly energy conserving. Its conservation and
convergence properties for these parameters remain unverified in this phase.
The nonzero initial COM velocity is retained; no recentering or numerical correction
has been applied. Long-term relaxation, close-encounter error, and softening sensitivity
remain open. No tests normalize, clip, or otherwise correct evolving quantities.

The legacy Jeans estimate is origin-dependent, assigns all mass to a 68th-percentile
radius, and resets that radius to the initial radius below softening. Colors/sizes
are radial illustration, not local density or temperature. Stage labels are elapsed
time thresholds, and the displayed percentage is 100 t/t_ff, not a contraction
measurement. README instability inequalities are reversed. Those historical behaviors
are documented, not corrected, in Phase 0.

Phase 1 should verify analytic forces and two-body orbits, no-gravity motion,
softened potential energy, momentum and angular momentum, COM evolution, timestep
convergence, softening sensitivity, and short collapse benchmarks with predefined
tolerances. Phase 1 requires explicit authorization. No Phase 1 implementation has
been started.

## Files changed

- Updated: `README.md` (status notice and safe test command only).
- Added: `.gitignore`, `requirements-baseline.txt`, `baseline/manifest.json`,
  `baseline/reference.json`, `tests/baseline_support.py`, `tests/test_baseline.py`,
  and this document.
- Removed: none.
- Preserved byte for byte: `starformation.py`, `star_formation.mp4`, `star_formation.gif`.
