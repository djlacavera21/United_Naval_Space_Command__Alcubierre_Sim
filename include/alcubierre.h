#ifndef ALCUBIERRE_H
#define ALCUBIERRE_H

#include <stddef.h>

#define ALC_VERSION "0.2.0"
#define ALC_C 299792458.0
#define ALC_G 6.67430e-11
#define ALC_PI 3.14159265358979323846

typedef struct { double x, y, z; } AlcVec3;

typedef struct {
    double radius_m, sigma_per_m, warp_factor;
    AlcVec3 direction;
    double ramp_up_s, cruise_s, ramp_down_s;
    double duration_s, sample_dt_s;
    unsigned energy_intervals;
} AlcConfig;

typedef struct {
    double time_s, throttle, speed_m_s, distance_m;
    AlcVec3 position_m;
} AlcState;

typedef struct {
    double shape;
    AlcVec3 gradient_per_m;
    double energy_density_j_m3;
    double expansion_per_s;
    /* Coordinates (t, x, y, z), with t in seconds. */
    double metric[4][4];
} AlcField;

typedef struct {
    double energy_j;
    double relative_difference;
    double radial_min_m, radial_max_m;
    unsigned intervals;
    int converged;
} AlcEnergy;

AlcConfig alc_default_config(void);
/* Returns 1 on success, 0 with a readable diagnostic otherwise. */
int alc_validate(const AlcConfig *config, char *error, size_t error_size);
AlcVec3 alc_unit(AlcVec3 vector);

/* Remaining functions require a successfully validated configuration.
 * Times/offsets must be finite; radii must be finite and nonnegative.
 * Field throttle must be in [0, 1]; sample indices must be below sample_count.
 * Fields take offsets relative to the bubble, avoiding subtraction of
 * astronomical coordinates when resolving a small bubble wall. */
double alc_shape(const AlcConfig *config, double radius_m);
double alc_shape_derivative(const AlcConfig *config, double radius_m);
double alc_throttle(const AlcConfig *config, double time_s);
double alc_profile_integral(const AlcConfig *config, double time_s);
AlcState alc_state_at(const AlcConfig *config, double time_s);
AlcField alc_field_at(const AlcConfig *config, double throttle, AlcVec3 offset_m);
AlcEnergy alc_energy_at(const AlcConfig *config, double throttle);
/* Includes t=0 and the requested endpoint, each exactly once. */
size_t alc_sample_count(const AlcConfig *config);
double alc_sample_time(const AlcConfig *config, size_t index);

#endif
