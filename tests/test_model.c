#include "alcubierre.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static unsigned checks = 0;

#define CHECK(condition) do { \
    ++checks; \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: failed: %s\n", __FILE__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
} while (0)

static int close_to(double actual, double expected, double relative, double absolute)
{
    return isfinite(actual) &&
        fabs(actual - expected) <= absolute + relative * fabs(expected);
}

static void test_shape(void)
{
    AlcConfig c = alc_default_config();
    CHECK(alc_shape(&c, 0.0) == 1.0);
    CHECK(alc_shape_derivative(&c, 0.0) == 0.0);
    double previous = 1.0;
    for (unsigned i = 0; i <= 500; ++i) {
        double r = 0.5 * i;
        double f = alc_shape(&c, r);
        CHECK(f >= 0.0 && f <= previous + 1e-15);
        CHECK(alc_shape_derivative(&c, r) <= 0.0);
        previous = f;
    }
    /* Compare the stable formulation to the original normalized tanh expression. */
    for (unsigned i = 1; i <= 25; ++i) {
        double r = 5.0 * i;
        double direct = (tanh(c.sigma_per_m * (r + c.radius_m)) -
                         tanh(c.sigma_per_m * (r - c.radius_m))) /
                         (2.0 * tanh(c.sigma_per_m * c.radius_m));
        CHECK(close_to(alc_shape(&c, r), direct, 1e-12, 1e-15));
    }
    for (unsigned i = 0; i < 7; ++i) {
        double r = 70.0 + 10.0 * i, h = 1e-3;
        double finite_difference = (alc_shape(&c, r + h) - alc_shape(&c, r - h)) / (2.0 * h);
        CHECK(close_to(alc_shape_derivative(&c, r), finite_difference, 2e-7, 1e-11));
    }
    c.sigma_per_m = 1e-8; /* sigma*R = 1e-6: broad-wall limiting profile */
    CHECK(close_to(alc_shape(&c, 1.0 / c.sigma_per_m),
                   1.0 / (cosh(1.0) * cosh(1.0)), 1e-10, 0.0));
    c.sigma_per_m = 1e4; /* sigma*R = 1e6: very thin wall */
    CHECK(alc_shape(&c, 0.0) == 1.0);
    CHECK(close_to(alc_shape(&c, c.radius_m), 0.5, 1e-14, 0.0));
    CHECK(close_to(alc_shape_derivative(&c, c.radius_m), -0.5 * c.sigma_per_m, 1e-14, 0.0));
    CHECK(alc_shape(&c, DBL_MAX) == 0.0);
    CHECK(alc_shape_derivative(&c, DBL_MAX) == 0.0);
}

static void test_fields(void)
{
    AlcConfig c = alc_default_config();
    AlcVec3 center = {0.0, 0.0, 0.0};
    AlcField flat = alc_field_at(&c, 0.0, center);
    CHECK(flat.metric[0][0] == -ALC_C * ALC_C);
    CHECK(flat.energy_density_j_m3 == 0.0);
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 4; ++j)
            CHECK(flat.metric[i][j] == (i == j ? (i == 0 ? -ALC_C * ALC_C : 1.0) : 0.0));

    AlcField axial = alc_field_at(&c, 1.0, (AlcVec3){100.0, 0.0, 0.0});
    AlcField transverse = alc_field_at(&c, 1.0, (AlcVec3){0.0, 100.0, 0.0});
    AlcField behind = alc_field_at(&c, 1.0, (AlcVec3){-100.0, 0.0, 0.0});
    CHECK(axial.energy_density_j_m3 == 0.0);
    CHECK(transverse.energy_density_j_m3 < 0.0);
    CHECK(axial.expansion_per_s < 0.0 && behind.expansion_per_s > 0.0);
    CHECK(transverse.expansion_per_s == 0.0);
    double original_rho = transverse.energy_density_j_m3;
    c.warp_factor *= 2.0;
    AlcField faster = alc_field_at(&c, 1.0, (AlcVec3){0.0, 100.0, 0.0});
    CHECK(close_to(faster.energy_density_j_m3, 4.0 * original_rho, 1e-14, 0.0));
    c = alc_default_config();
    c.direction = (AlcVec3){0.0, 9.0, 0.0};
    AlcField rotated = alc_field_at(&c, 1.0, (AlcVec3){100.0, 0.0, 0.0});
    CHECK(close_to(rotated.energy_density_j_m3, original_rho, 1e-14, 0.0));

    c.direction = (AlcVec3){1.0, 2.0, 3.0};
    AlcVec3 unit = alc_unit(c.direction);
    AlcField field = alc_field_at(&c, 0.75, (AlcVec3){55.0, 70.0, 11.0});
    double n[3] = {unit.x, unit.y, unit.z};
    double shift_squared = 0.0;
    for (size_t i = 1; i < 4; ++i)
        shift_squared += field.metric[0][i] * field.metric[0][i];
    /* Schur complement = determinant for the identity spatial block. */
    CHECK(close_to(field.metric[0][0] - shift_squared, -ALC_C * ALC_C, 1e-14, 0.0));
    for (int sign = -1; sign <= 1; sign += 2) {
        double ds = field.metric[0][0] / (ALC_C * ALC_C);
        for (size_t i = 0; i < 3; ++i) {
            double velocity_c = (c.warp_factor * 0.75 * field.shape + sign) * n[i];
            ds += 2.0 * field.metric[0][i + 1] / ALC_C * velocity_c + velocity_c * velocity_c;
        }
        CHECK(fabs(ds) < 1e-13); /* Coordinate light-cone branches are null. */
    }
    AlcField at_center = alc_field_at(&c, 1.0, center);
    double ds = at_center.metric[0][0] / (ALC_C * ALC_C);
    for (size_t i = 0; i < 3; ++i) {
        double velocity_c = c.warp_factor * n[i];
        ds += 2.0 * at_center.metric[0][i + 1] / ALC_C * velocity_c + velocity_c * velocity_c;
    }
    CHECK(close_to(ds, -1.0, 1e-13, 0.0));
}

static void test_motion(void)
{
    AlcConfig c = alc_default_config();
    CHECK(alc_throttle(&c, 0.0) == 0.0);
    CHECK(alc_throttle(&c, 0.125) == 0.5);
    CHECK(alc_throttle(&c, 0.25) == 1.0);
    CHECK(alc_throttle(&c, 1.25) == 1.0);
    CHECK(alc_throttle(&c, 1.375) == 0.5);
    CHECK(alc_throttle(&c, 1.5) == 0.0);
    CHECK(alc_profile_integral(&c, 0.0) == 0.0);
    CHECK(close_to(alc_profile_integral(&c, 1.6), 1.25, 1e-15, 0.0));
    AlcState final = alc_state_at(&c, 1.6);
    CHECK(close_to(final.position_m.x, 749481145.0, 1e-14, 0.0));
    CHECK(final.speed_m_s == 0.0);
    CHECK(final.position_m.y == 0.0 && final.position_m.z == 0.0);
    for (unsigned i = 1; i < 160; ++i) {
        double t = 0.01 * i, h = 1e-5;
        double numerical_velocity = (alc_state_at(&c, t + h).distance_m -
                                     alc_state_at(&c, t - h).distance_m) / (2.0 * h);
        CHECK(close_to(numerical_velocity, alc_state_at(&c, t).speed_m_s, 1e-6, 1.0));
    }
    c.sample_dt_s = 0.137; /* Sampling must not alter the exact trajectory. */
    CHECK(alc_state_at(&c, 1.6).position_m.x == final.position_m.x);
    c.direction = (AlcVec3){0.0, 0.0, 5.0};
    c.warp_factor = -2.0;
    final = alc_state_at(&c, 1.6);
    CHECK(final.position_m.x == 0.0 && final.position_m.y == 0.0);
    CHECK(close_to(final.position_m.z, -749481145.0, 1e-14, 0.0));
    c.cruise_s = 0.0;
    CHECK(close_to(alc_profile_integral(&c, 5.0), 0.25, 1e-14, 0.0));
}

static void test_energy(void)
{
    AlcConfig c = alc_default_config();
    AlcEnergy e = alc_energy_at(&c, 1.0);
    CHECK(e.converged && e.energy_j < 0.0);
    CHECK(e.relative_difference <= 1e-6);
    CHECK(close_to(alc_energy_at(&c, 0.5).energy_j, 0.25 * e.energy_j, 1e-14, 0.0));
    CHECK(alc_energy_at(&c, 0.0).energy_j == 0.0);
    c.direction = (AlcVec3){1.0, 7.0, 9.0};
    CHECK(alc_energy_at(&c, 1.0).energy_j == e.energy_j);
    c.radius_m *= 2.0;
    c.sigma_per_m /= 2.0;
    CHECK(close_to(alc_energy_at(&c, 1.0).energy_j, 2.0 * e.energy_j, 1e-13, 0.0));

    c = alc_default_config();
    c.sigma_per_m = 10.0; /* Independent thin-wall asymptotic integral. */
    e = alc_energy_at(&c, 1.0);
    double thin_wall = -(ALC_C * ALC_C) * (ALC_C * ALC_C) / (36.0 * ALC_G) *
        c.warp_factor * c.warp_factor * c.radius_m * c.radius_m * c.sigma_per_m;
    CHECK(e.converged);
    CHECK(close_to(e.energy_j, thin_wall, 1e-5, 0.0));
    c.sigma_per_m = 1e-8;
    CHECK(alc_energy_at(&c, 1.0).converged);
    c.sigma_per_m = 1e4;
    CHECK(alc_energy_at(&c, 1.0).converged);
}

static void test_validation_and_sampling(void)
{
    AlcConfig c = alc_default_config();
    char error[200];
    CHECK(alc_validate(&c, error, sizeof(error)));
    CHECK(alc_sample_count(&c) == 1601);
    c.duration_s = 0.3;
    c.sample_dt_s = 0.1;
    CHECK(alc_sample_count(&c) == 4);
    CHECK(alc_sample_time(&c, 0) == 0.0);
    CHECK(alc_sample_time(&c, 3) == 0.3);
    c.sample_dt_s = 0.07;
    CHECK(alc_sample_count(&c) == 6);
    for (size_t i = 1; i < alc_sample_count(&c); ++i)
        CHECK(alc_sample_time(&c, i) > alc_sample_time(&c, i - 1));
    CHECK(alc_sample_time(&c, 5) == 0.3);
    c.duration_s = 0.0;
    CHECK(alc_sample_count(&c) == 1 && alc_sample_time(&c, 0) == 0.0);
    c.duration_s = 1e-100;
    c.sample_dt_s = 1e300;
    CHECK(alc_sample_count(&c) == 2);
    c = alc_default_config();
    c.direction = (AlcVec3){DBL_MAX, DBL_MAX, 0.0};
    CHECK(alc_validate(&c, error, sizeof(error)));
    AlcVec3 n = alc_unit(c.direction);
    CHECK(close_to(hypot(n.x, n.y), 1.0, 1e-14, 0.0));
    c.direction = (AlcVec3){DBL_MIN, 0.0, 0.0};
    CHECK(alc_unit(c.direction).x == 1.0);

    double *fields[] = {&c.radius_m, &c.sigma_per_m, &c.warp_factor,
        &c.direction.x, &c.ramp_up_s, &c.cruise_s, &c.ramp_down_s,
        &c.duration_s, &c.sample_dt_s};
    for (size_t i = 0; i < sizeof(fields) / sizeof(fields[0]); ++i) {
        double old = *fields[i];
        *fields[i] = NAN;
        CHECK(!alc_validate(&c, error, sizeof(error)) && error[0] != '\0');
        *fields[i] = old;
    }
    c.direction = (AlcVec3){0.0, 0.0, 0.0};
    CHECK(!alc_validate(&c, error, sizeof(error)));
    c = alc_default_config();
    c.sample_dt_s = 0.0;
    CHECK(!alc_validate(&c, error, sizeof(error)));
    c = alc_default_config();
    c.energy_intervals = 33;
    CHECK(!alc_validate(&c, error, sizeof(error)));
    c = alc_default_config();
    c.ramp_up_s = 0.0;
    CHECK(!alc_validate(&c, error, sizeof(error)));
    CHECK(!alc_validate(NULL, NULL, 0));
}

int main(void)
{
    test_shape();
    test_fields();
    test_motion();
    test_energy();
    test_validation_and_sampling();
    printf("Model: %u numerical assertions passed\n", checks);
    return EXIT_SUCCESS;
}
