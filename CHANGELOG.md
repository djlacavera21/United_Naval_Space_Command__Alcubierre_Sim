# Changelog

## 0.2.0 — 2026-09-14

### Added

- Buildable C11 model library and command-line simulator.
- Make targets for compilation, model tests, and CLI regression tests.
- Runtime parameters for geometry, direction, speed, timing, and quadrature.
- JSON run summaries and optional comoving 2D field CSV exports.
- Analytic shape gradients, expansion diagnostics, and radial energy convergence.
- GCC/Clang/macOS CI and address/undefined-behavior sanitizer checks.
- Equations, units, numerical limits, and documented next milestones.

### Corrected

- Telemetry previously advanced the clock and position before writing a throttle
  calculated at the prior timestamp. Every value now belongs to the stated time.
- Exact throttle integration replaces timestep-dependent Euler updates for the
  prescribed center trajectory.
- Stable exponential profile evaluation preserves broad/thin wall behavior.
- Energy density now uses the gradient perpendicular to the travel direction,
  replacing the original full-gradient toy proxy.
- Radial integration resolves the wall and evaluates the angular integral
  analytically, replacing repeated finite-difference angular sampling.
- Existing output files are preserved; failed writes/open operations return
  nonzero and discard new files created during the unsuccessful run.

### Compatibility

- The original alcubierre_sim.txt is retained as an unbuilt historical prototype.
- Build sources now live under src/ and include/.
- CSV schema now uses explicit units and samples the transverse wall.
- Energy results intentionally differ from the original toy estimate.
