#include "alcubierre.h"

#include <float.h>
#include <math.h>
#include <stdio.h>

static double length(AlcVec3 v)
{
    return hypot(hypot(v.x, v.y), v.z);
}

static double dot(AlcVec3 a, AlcVec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static AlcVec3 multiply(AlcVec3 v, double factor)
{
    AlcVec3 result = {v.x * factor, v.y * factor, v.z * factor};
    return result;
}

static AlcVec3 cross(AlcVec3 a, AlcVec3 b)
{
    AlcVec3 result = {a.y * b.z - a.z * b.y,
                      a.z * b.x - a.x * b.z,
                      a.x * b.y - a.y * b.x};
    return result;
}

AlcVec3 alc_unit(AlcVec3 v)
{
    double largest = fmax(fabs(v.x), fmax(fabs(v.y), fabs(v.z)));
    if (!(largest > 0.0) || !isfinite(largest)) {
        AlcVec3 zero = {0.0, 0.0, 0.0};
        return zero;
    }
    /* Divide components individually: 1/largest could overflow. */
    AlcVec3 scaled = {v.x / largest, v.y / largest, v.z / largest};
    return multiply(scaled, 1.0 / length(scaled));
}

AlcConfig alc_default_config(void)
{
    AlcConfig c = {100.0, 0.10, 2.0, {1.0, 0.0, 0.0},
                   0.25, 1.0, 0.25, 1.60, 0.001, 1024};
    return c;
}

static int invalid(char *error, size_t capacity, const char *message)
{
    if (error != NULL && capacity > 0)
        (void)snprintf(error, capacity, "%s", message);
    return 0;
}

int alc_validate(const AlcConfig *c, char *error, size_t capacity)
{
    if (c == NULL)
        return invalid(error, capacity, "configuration is required");
    if (!isfinite(c->radius_m) || c->radius_m < 1e-6 || c->radius_m > 1e9)
        return invalid(error, capacity, "radius must be in [1e-6, 1e9] meters");
    if (!isfinite(c->sigma_per_m) || c->sigma_per_m <= 0.0)
        return invalid(error, capacity, "sigma must be finite and positive");
    double a = c->sigma_per_m * c->radius_m;
    if (!isfinite(a) || a < 1e-6 || a > 1e6)
        return invalid(error, capacity, "sigma * radius must be in [1e-6, 1e6]");
    if (!isfinite(c->warp_factor) || fabs(c->warp_factor) > 1000.0)
        return invalid(error, capacity, "warp factor must be in [-1000, 1000]");
    if (!isfinite(c->direction.x) || !isfinite(c->direction.y) ||
        !isfinite(c->direction.z) ||
        (c->direction.x == 0.0 && c->direction.y == 0.0 && c->direction.z == 0.0))
        return invalid(error, capacity, "direction must be finite and nonzero");
    if (!isfinite(c->ramp_up_s) || c->ramp_up_s < 1e-9 ||
        !isfinite(c->ramp_down_s) || c->ramp_down_s < 1e-9 ||
        !isfinite(c->cruise_s) || c->cruise_s < 0.0 ||
        c->ramp_up_s + c->cruise_s + c->ramp_down_s > 1e9)
        return invalid(error, capacity, "ramps must be >= 1e-9 s, cruise >= 0, profile <= 1e9 s");
    if (!isfinite(c->duration_s) || c->duration_s < 0.0 || c->duration_s > 1e9)
        return invalid(error, capacity, "duration must be in [0, 1e9] seconds");
    if (!isfinite(c->sample_dt_s) || c->sample_dt_s <= 0.0 ||
        c->duration_s / c->sample_dt_s > 1000000.0)
        return invalid(error, capacity, "dt must be positive and request at most 1000000 intervals");
    if (c->energy_intervals < 32 || c->energy_intervals > 65536 ||
        c->energy_intervals % 2 != 0)
        return invalid(error, capacity, "energy intervals must be even and in [32, 65536]");
    if (error != NULL && capacity > 0)
        error[0] = '\0';
    return 1;
}

/* Algebraically equivalent to the normalized tanh profile. All exponentials
 * have nonpositive arguments, avoiding overflow and far-field cancellation. */
double alc_shape(const AlcConfig *c, double r)
{
    double a = c->sigma_per_m * c->radius_m;
    double b = c->sigma_per_m * fabs(r);
    if (!isfinite(b))
        return 0.0;
    double q = exp(-2.0 * a);
    double numerator = exp(-2.0 * fmax(b - a, 0.0)) * (1.0 + q) * (1.0 + q);
    double denominator = (1.0 + exp(-2.0 * (a + b))) *
                         (1.0 + exp(-2.0 * fabs(a - b)));
    return fmin(1.0, numerator / denominator);
}

double alc_shape_derivative(const AlcConfig *c, double r)
{
    double a = c->sigma_per_m * c->radius_m;
    double b = c->sigma_per_m * fabs(r);
    if (b == 0.0 || !isfinite(b))
        return 0.0;
    double m = fmax(a, b);
    double eb = exp(2.0 * (b - m));
    double ea = exp(2.0 * (a - m));
    /* Stable form of tanh(a+b) + tanh(b-a), including near r=0. */
    double slope = 2.0 * eb * (-expm1(-4.0 * b)) /
        (eb * (1.0 + exp(-4.0 * b)) + ea * (1.0 + exp(-4.0 * a)));
    return -c->sigma_per_m * alc_shape(c, r) * slope;
}

static double smoothstep(double x)
{
    return x * x * (3.0 - 2.0 * x);
}

static double smoothstep_integral(double x)
{
    return x * x * x * (1.0 - 0.5 * x);
}

double alc_throttle(const AlcConfig *c, double t)
{
    if (t <= 0.0)
        return 0.0;
    if (t < c->ramp_up_s)
        return smoothstep(t / c->ramp_up_s);
    t -= c->ramp_up_s;
    if (t < c->cruise_s)
        return 1.0;
    t -= c->cruise_s;
    if (t < c->ramp_down_s)
        return 1.0 - smoothstep(t / c->ramp_down_s);
    return 0.0;
}

double alc_profile_integral(const AlcConfig *c, double t)
{
    if (t <= 0.0)
        return 0.0;
    if (t < c->ramp_up_s)
        return c->ramp_up_s * smoothstep_integral(t / c->ramp_up_s);
    double total = 0.5 * c->ramp_up_s;
    t -= c->ramp_up_s;
    if (t < c->cruise_s)
        return total + t;
    total += c->cruise_s;
    t -= c->cruise_s;
    if (t < c->ramp_down_s) {
        double x = t / c->ramp_down_s;
        return total + c->ramp_down_s * (x - smoothstep_integral(x));
    }
    return total + 0.5 * c->ramp_down_s;
}

AlcState alc_state_at(const AlcConfig *c, double t)
{
    AlcState state;
    state.time_s = t;
    state.throttle = alc_throttle(c, t);
    state.speed_m_s = ALC_C * c->warp_factor * state.throttle;
    state.distance_m = ALC_C * c->warp_factor * alc_profile_integral(c, t);
    state.position_m = multiply(alc_unit(c->direction), state.distance_m);
    return state;
}

AlcField alc_field_at(const AlcConfig *c, double throttle, AlcVec3 offset)
{
    AlcField field = {0};
    double r = length(offset);
    AlcVec3 direction = alc_unit(c->direction);
    double derivative = alc_shape_derivative(c, r);
    field.shape = alc_shape(c, r);
    if (r > 0.0) {
        AlcVec3 radial = {offset.x / r, offset.y / r, offset.z / r};
        field.gradient_per_m = multiply(radial, derivative);
    }
    /* Eq. (19): only the gradient perpendicular to the travel direction
     * contributes. A cross product avoids subtracting nearly equal squares. */
    AlcVec3 transverse = cross(direction, field.gradient_per_m);
    double beta = c->warp_factor * throttle;
    double coefficient = (ALC_C * ALC_C) * (ALC_C * ALC_C) /
                         (32.0 * ALC_PI * ALC_G);
    field.energy_density_j_m3 = -coefficient * beta * beta * dot(transverse, transverse);
    double speed = ALC_C * beta;
    field.expansion_per_s = speed * dot(direction, field.gradient_per_m);
    double n[3] = {direction.x, direction.y, direction.z};
    field.metric[0][0] = -ALC_C * ALC_C + speed * speed * field.shape * field.shape;
    for (size_t i = 0; i < 3; ++i) {
        field.metric[0][i + 1] = -speed * field.shape * n[i];
        field.metric[i + 1][0] = field.metric[0][i + 1];
        field.metric[i + 1][i + 1] = 1.0;
    }
    return field;
}

/* Integrate in u=sigma*(r-R) so a thin wall is always resolved.
 * The angular integral of sin(theta)^2 is exactly 8*pi/3. */
static double radial_integral(const AlcConfig *c, unsigned intervals)
{
    double a = c->sigma_per_m * c->radius_m;
    double lower = fmax(-a, -12.0);
    double du = (12.0 - lower) / intervals;
    double sum = 0.0;
    for (unsigned i = 0; i <= intervals; ++i) {
        double u = lower + du * i;
        double r = fmax(0.0, (a + u) / c->sigma_per_m);
        double derivative = alc_shape_derivative(c, r);
        double value = (r * derivative) * (r * derivative) / c->sigma_per_m;
        double weight = (i == 0 || i == intervals) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
        sum += weight * value;
    }
    return sum * du / 3.0;
}

AlcEnergy alc_energy_at(const AlcConfig *c, double throttle)
{
    AlcEnergy result = {0};
    result.radial_min_m = fmax(0.0, c->radius_m - 12.0 / c->sigma_per_m);
    result.radial_max_m = c->radius_m + 12.0 / c->sigma_per_m;
    result.intervals = 2 * c->energy_intervals;
    double beta = c->warp_factor * throttle;
    if (beta == 0.0) {
        result.converged = 1;
        return result;
    }
    double coarse = radial_integral(c, c->energy_intervals);
    double fine = radial_integral(c, result.intervals);
    double coefficient = (ALC_C * ALC_C) * (ALC_C * ALC_C) / (12.0 * ALC_G);
    result.energy_j = -coefficient * beta * beta * fine;
    result.relative_difference = fabs(fine - coarse) / fmax(fabs(fine), DBL_MIN);
    result.converged = isfinite(result.energy_j) &&
        isfinite(result.relative_difference) && result.relative_difference <= 1e-6;
    return result;
}

size_t alc_sample_count(const AlcConfig *c)
{
    if (c->duration_s == 0.0)
        return 1;
    size_t intervals = (size_t)ceil(c->duration_s / c->sample_dt_s);
    if (intervals == 0)
        intervals = 1;
    if (intervals > 1 &&
        fabs((intervals - 1) * c->sample_dt_s - c->duration_s) <=
        8.0 * DBL_EPSILON * c->duration_s)
        --intervals;
    return intervals + 1;
}

double alc_sample_time(const AlcConfig *c, size_t index)
{
    return index + 1 == alc_sample_count(c) ? c->duration_s : index * c->sample_dt_s;
}
