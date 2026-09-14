# Model and numerical conventions

## Prescribed spacetime

The model follows Miguel Alcubierre, *The warp drive: hyper-fast travel within
general relativity*, Classical and Quantum Gravity 11 (1994), L73–L77.
[Original paper](https://arxiv.org/abs/gr-qc/0009013), equations (6), (8), and (19).

We use SI units and coordinates (t,x,y,z), with a fixed unit direction n and
signed bubble speed v(t). The prescribed line element is

$$
ds^2=-c^2dt^2+\left|d\mathbf{x}-v(t)f(r)\mathbf{n}\,dt\right|^2,
\qquad r=|\mathbf{x}-\mathbf{x}_s(t)|.
$$

The shift is -vf n. Therefore

$$
g_{tt}=-c^2+v^2f^2,\qquad g_{ti}=-vf n_i,\qquad g_{ij}=\delta_{ij}.
$$

The code reports these components directly. The temporal coordinate is seconds,
not ct; the temporal components consequently carry units. At zero speed the
metric reduces to Minkowski spacetime. At the bubble center f=1; substituting
the prescribed center velocity into the metric yields ds²/dt²=-c².
These relationships are numerical regression checks, not demonstrations of a
propulsion mechanism.

## Shape and derivatives

The normalized profile is

$$
f(r)=\frac{\tanh[\sigma(r+R)]-\tanh[\sigma(r-R)]}{2\tanh(\sigma R)}.
$$

Let a=σR and b=σr. We evaluate the equivalent expression

$$
f(r)=
\frac{e^{-2\max(b-a,0)}(1+e^{-2a})^2}
{(1+e^{-2(a+b)})(1+e^{-2|a-b|})}.
$$

This avoids subtraction of saturated tanh values and overflow in hyperbolic
functions. The analytic derivative is evaluated with a stable exponential form
of the factor in brackets:

$$
f'(r)=-\sigma f(r)\,[\tanh(a+b)+\tanh(b-a)],
\qquad
\nabla f=f'(r)\frac{\mathbf{r}}{r}.
$$

The origin has exactly zero gradient. C API field coordinates are already
relative to the bubble center, avoiding cancellation when a small bubble moves
a large distance. The API requires validated configurations and finite arguments;
radial distances and sample indices must satisfy the header's preconditions.

## Energy and expansion

The Eulerian energy density for this prescribed metric, restoring SI units in
equation (19) of the paper, is

$$
\rho_E=-\frac{c^4}{32\pi G}\left(\frac{v}{c}\right)^2
\left|\mathbf{n}\times\nabla f\right|^2.
$$

Its units are J/m³. The original uploaded prototype used the entire squared
gradient as a qualitative proxy. Version 0.2.0 uses the transverse gradient, so
the energy density vanishes along the direction axis. As a result, the old
wall-energy column and old integrated-energy values are not directly compatible
with the new output.

For a spherical profile, angular integration gives

$$
E=-\frac{c^4}{12G}\left(\frac{v}{c}\right)^2
\int_0^\infty r^2[f'(r)]^2\,dr.
$$

This is an integral on the chosen flat spatial slice. It is not a conserved
global energy, hardware input power, fuel consumption, or a self-consistent
matter model. Dividing by c² reports a signed mass equivalent only. The
expansion diagnostic is the spatial divergence of the Eulerian velocity:

$$
\Theta=v\,\mathbf{n}\cdot\nabla f.
$$

For positive speed, it is negative in front and positive behind. The code does
not calculate the complete stress-energy tensor, tidal forces, horizons,
quantum inequalities, or the evolution of a matter source.

## Quadrature

The implementation changes variables to u=σ(r-R) and integrates from
max(-σR,-12) to +12. This confines the mesh to the wall even when σR is large.
The radial bounds are recorded in the JSON summary. The omitted inner and outer
tails are exponentially suppressed; no rigorous tail bound is reported.

Composite Simpson quadrature runs with N and 2N even intervals. The finer
result is returned. The convergence diagnostic is their relative difference:

$$
\delta=\frac{|I_{2N}-I_N|}{\max(|I_{2N}|,\mathrm{DBL\_MIN})}.
$$

The CLI requires δ <= 1e-6 before creating output files. Increasing
--energy-intervals can resolve a failure. This is an empirical mesh comparison,
not a rigorous error bound. Floating-point and truncation error are separate.

Full-throttle quadrature is performed once per run. At other times its value is
multiplied by throttle squared. The zero-speed case returns exact zero without
quadrature. Tests check the thin-wall asymptotic result

$$
E\sim-\frac{c^4}{36G}\left(\frac{v}{c}\right)^2R^2\sigma
\quad(\sigma R\gg1).
$$

## Motion and sampling

The schedule has cubic smoothstep S(q)=3q²-2q³ on ramp-up, a constant cruise,
and 1-S(q) on ramp-down. Its primitive F(q)=q³-q⁴/2 allows direct evaluation of

$$
\mathbf{x}_s(t)=c\,w\,\mathbf{n}\int_0^t T(s)\,ds
$$

for requested signed warp factor w and throttle T. After the complete schedule,

$$
\int T(s)\,ds=\tfrac12 t_{\mathrm{up}}+t_{\mathrm{cruise}}+
\tfrac12 t_{\mathrm{down}}.
$$

The timestep controls output sampling only. The program exports t=0 and the
requested endpoint exactly once, snapping roundoff-scale endpoint differences.
Each row evaluates throttle, position, velocity, and fields at the same time.
No numerical trajectory solver is claimed.

The program currently fixes the direction, radius, and steepness throughout a
run. Changing those in time would require extending the model and its tests.

## Verification

Run **make test** to compile and execute the C assertions and Python tests of
the compiled CLI. The tests cover analytic limits, finite-difference derivative
cross-checks, metric light-cone consistency, rotational/scaling relations,
convergence, exact motion, CLI rejection paths, CSV/JSON output, and preservation
of existing files. CI runs GCC and Clang on Linux, Clang on macOS, and Clang with
AddressSanitizer and UndefinedBehaviorSanitizer.
