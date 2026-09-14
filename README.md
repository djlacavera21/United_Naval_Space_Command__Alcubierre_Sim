# United Naval Space Command — Alcubierre Simulator

Version **0.2.0** turns the original C prototype into a configurable C11 numerical
simulator with a reusable model API, reproducible telemetry, and field diagnostics.

This project evaluates a **prescribed Alcubierre spacetime metric**. It does not
operate propulsion hardware or demonstrate a physically realizable warp drive.
The motion is a specified bubble-center trajectory; no general geodesic or
Einstein-equation evolution solver is included.

## Build and run

Requires a C11 compiler, Make, and the system math library. Tests additionally
require Python 3.9 or later and use only its standard library.

~~~sh
make
./build/alcubierre_sim
make test
~~~

The default run writes 1,601 samples from 0 through 1.6 seconds. It ramps to a
coordinate speed of 2c, cruises, and ramps down. Its final coordinate displacement
is **749,481,145 meters**.

Output files are created exclusively: existing files are preserved and the run
fails with a diagnostic. Choose new paths or remove previous outputs yourself
before rerunning. If opening or writing an output fails, newly created outputs
from that run are removed.

A custom run with a field snapshot:

~~~sh
./build/alcubierre_sim \
  --radius 100 --sigma 0.1 --warp 2 \
  --direction 1 0 0 \
  --ramp-up 0.25 --cruise 1 --ramp-down 0.25 \
  --duration 1.6 --dt 0.005 \
  --output experiment.csv --summary experiment.json \
  --slice experiment_field.csv --slice-time 0.75 \
  --slice-points 101 --slice-extent 160
~~~

Use `./build/alcubierre_sim --help` for the complete option list. Negative warp
factors reverse motion relative to the configured direction. A zero warp factor
provides the flat-spacetime control case. A duration of zero exports just the
initial state; durations shorter than the throttle schedule export a partial run.

## What is implemented

- Stable bubble shape evaluation and analytic radial derivatives.
- Metric components for any fixed 3D travel direction, in coordinates (t,x,y,z).
- Eulerian expansion and directional energy density for the prescribed metric.
- Full angular reduction of the energy integral, with Simpson quadrature in a
  wall-relative coordinate and a coarse/fine convergence diagnostic.
- Exact integration of the cubic throttle profile, independent of telemetry dt.
- Matching timestamps, velocities, positions, and fields, including both endpoints.
- CSV trajectory output, JSON configuration/results, and optional 2D field slices.
- Finite-input validation, bounded sampling, checked output writes, and regression
  tests under GCC, Clang, macOS, and address/undefined-behavior sanitizers in CI.

## Outputs

| File | Contents |
| --- | --- |
| `warp_telemetry.csv` | Time, throttle, signed coordinate speed, position, transverse-wall energy density, and instantaneous integrated energy |
| `warp_summary.json` | Schema/version, normalized configuration, endpoint, full-throttle energy, quadrature convergence, and slice settings |
| Optional `--slice FILE` | Comoving longitudinal/transverse coordinates, shape, energy density, expansion, temporal metric component, and parallel shift |

All dimensional CSV columns include units in their names. Telemetry and JSON use
17 significant digits. Energy is in joules, density in joules per cubic meter,
and the E/c² diagnostic is a signed mass equivalent in kilograms.

Slice coordinates lie in a meridional plane containing the configured travel
axis. They are offsets from the bubble center, preserving wall resolution even
after the bubble has traveled a large coordinate distance. Positive longitudinal
coordinates follow `--direction`; a negative warp factor reverses which side is
the front. The CSV is a sampled field, so a coarse slice may miss a thin wall;
adjust the extent and point count when inspecting steep profiles.

The reported full-throttle energy is evaluated at the requested warp factor even
if a partial run never reaches full throttle. The instantaneous telemetry energy
scales with throttle squared. It is a spatial energy integral of this metric,
**not fuel consumption, engine power, or an engineering energy budget**.

## Numerical limits and interpretation

| Parameter | Accepted range |
| --- | --- |
| Radius | 1e-6 to 1e9 m |
| Shape steepness | Positive, with 1e-6 <= sigma * radius <= 1e6 |
| Signed warp factor | -1000 to 1000 |
| Direction | Finite and nonzero; normalized internally |
| Each ramp | At least 1e-9 s |
| Cruise | Nonnegative; total throttle schedule at most 1e9 s |
| Observation duration | 0 to 1e9 s |
| Sampling interval | Positive; at most 1,000,000 intervals |
| Coarse energy quadrature | Even interval count, 32 to 65536 |
| Slice | 3 to 501 points per axis; positive extent at most 1e17 m |

The numerical bounds are implementation limits, not physical feasibility limits.
Energy convergence compares N and 2N intervals and requires a relative difference
at most 1e-6. This diagnostic does not bound the exponentially suppressed
radial tails outside the integration domain, floating-point error, or model error.
Read [the equations and numerical conventions](docs/model.md) before interpreting
energy values.

## Source layout

- `include/alcubierre.h`: reusable C API and units.
- `src/alcubierre.c`: fields, motion, sampling, and energy quadrature.
- `src/main.c`: command-line interface and exports.
- `tests/`: model assertions and compiled-program regression tests.
- `alcubierre_sim.txt`: original prototype retained for provenance; not built.
- `.github/workflows/ci.yml`: compiler matrix and sanitizer checks.

To build without Make:

~~~sh
cc -O2 -std=c11 -Wall -Wextra -Wpedantic -Iinclude \
  src/main.c src/alcubierre.c -lm -o alcubierre_sim
~~~

## Next milestones

1. Null-ray integration with step-size control and constraint monitoring.
2. Independent tensor/curvature diagnostics against symbolic reference cases.
3. A field viewer and parameter-sweep tooling with convergence comparisons.

See [CHANGELOG.md](CHANGELOG.md) for changes from the uploaded prototype.
