#include "alcubierre.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    AlcConfig config;
    const char *output_path, *summary_path, *slice_path;
    double slice_time_s, slice_extent_m;
    unsigned slice_points;
    int slice_options;
} Options;

static void usage(void)
{
    puts("UNSC Alcubierre Simulator " ALC_VERSION "\n"
         "Numerical spacetime model; no physical propulsion control.\n\n"
         "Usage: alcubierre_sim [options]\n"
         "  --radius M             Bubble radius (default 100)\n"
         "  --sigma INV_M          Wall steepness (default 0.1)\n"
         "  --warp FACTOR          Signed coordinate speed / c (default 2)\n"
         "  --direction X Y Z      Nonzero travel vector (default 1 0 0)\n"
         "  --ramp-up S            Smooth acceleration duration (default 0.25)\n"
         "  --cruise S             Cruise duration (default 1)\n"
         "  --ramp-down S          Smooth deceleration duration (default 0.25)\n"
         "  --duration S           Observation duration (default 1.6)\n"
         "  --dt S                 Telemetry sample interval (default 0.001)\n"
         "  --energy-intervals N   Even coarse Simpson interval count (default 1024)\n"
         "  --output FILE          Telemetry CSV (default warp_telemetry.csv)\n"
         "  --summary FILE         Run JSON (default warp_summary.json)\n"
         "  --slice FILE           Optional comoving field CSV\n"
         "  --slice-time S         Snapshot time (default cruise midpoint, capped at duration)\n"
         "  --slice-points N       Points per dimension, 3..501 (default 101)\n"
         "  --slice-extent M       Half-width (default radius + 6/sigma)\n"
         "  --help                 Show usage\n"
         "  --version              Show version\n\n"
         "Output files must not already exist. Use new filenames for a new run.");
}

static int number(const char *text, double *value)
{
    char *end = NULL;
    errno = 0;
    double parsed = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !isfinite(parsed))
        return 0;
    *value = parsed;
    return 1;
}

static int integer(const char *text, unsigned *value)
{
    if (*text == '\0')
        return 0;
    unsigned result = 0;
    for (const char *p = text; *p != '\0'; ++p) {
        if (*p < '0' || *p > '9' || result > 100000U)
            return 0;
        result = result * 10U + (unsigned)(*p - '0');
    }
    *value = result;
    return 1;
}

/* 0=error, 1=run, 2=informational exit */
static int parse(int argc, char **argv, Options *o)
{
    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        if (strcmp(arg, "--help") == 0) {
            usage();
            return 2;
        }
        if (strcmp(arg, "--version") == 0) {
            puts(ALC_VERSION);
            return 2;
        }
        if (strcmp(arg, "--direction") == 0) {
            if (i + 3 >= argc ||
                !number(argv[i + 1], &o->config.direction.x) ||
                !number(argv[i + 2], &o->config.direction.y) ||
                !number(argv[i + 3], &o->config.direction.z)) {
                fputs("--direction requires three finite numbers\n", stderr);
                return 0;
            }
            i += 3;
            continue;
        }
        if (i + 1 >= argc) {
            fprintf(stderr, "missing value for %s\n", arg);
            return 0;
        }
        const char *value = argv[++i];
        if (strcmp(arg, "--output") == 0) o->output_path = value;
        else if (strcmp(arg, "--summary") == 0) o->summary_path = value;
        else if (strcmp(arg, "--slice") == 0) o->slice_path = value;
        else if (strcmp(arg, "--energy-intervals") == 0) {
            if (!integer(value, &o->config.energy_intervals)) {
                fputs("invalid --energy-intervals integer\n", stderr);
                return 0;
            }
        } else if (strcmp(arg, "--slice-points") == 0) {
            o->slice_options = 1;
            if (!integer(value, &o->slice_points)) {
                fputs("invalid --slice-points integer\n", stderr);
                return 0;
            }
        } else {
            double *target = NULL;
            if (strcmp(arg, "--radius") == 0) target = &o->config.radius_m;
            else if (strcmp(arg, "--sigma") == 0) target = &o->config.sigma_per_m;
            else if (strcmp(arg, "--warp") == 0) target = &o->config.warp_factor;
            else if (strcmp(arg, "--dt") == 0) target = &o->config.sample_dt_s;
            else if (strcmp(arg, "--duration") == 0) target = &o->config.duration_s;
            else if (strcmp(arg, "--ramp-up") == 0) target = &o->config.ramp_up_s;
            else if (strcmp(arg, "--cruise") == 0) target = &o->config.cruise_s;
            else if (strcmp(arg, "--ramp-down") == 0) target = &o->config.ramp_down_s;
            else if (strcmp(arg, "--slice-time") == 0) {
                target = &o->slice_time_s;
                o->slice_options = 1;
            } else if (strcmp(arg, "--slice-extent") == 0) {
                target = &o->slice_extent_m;
                o->slice_options = 1;
            }
            if (target == NULL) {
                fprintf(stderr, "unknown option: %s\n", arg);
                return 0;
            }
            if (!number(value, target)) {
                fprintf(stderr, "invalid finite number for %s: %s\n", arg, value);
                return 0;
            }
        }
    }
    return 1;
}

static AlcVec3 transverse_direction(AlcVec3 n)
{
    AlcVec3 vector;
    if (fabs(n.x) < 0.9) {
        vector = (AlcVec3){0.0, n.z, -n.y};
    } else {
        vector = (AlcVec3){-n.z, 0.0, n.x};
    }
    return alc_unit(vector);
}

static int telemetry(FILE *file, const AlcConfig *c, AlcEnergy peak)
{
    if (fprintf(file, "time_s,throttle,warp_factor,velocity_m_s,velocity_c,"
        "x_m,y_m,z_m,shape_center,transverse_wall_energy_density_j_m3,"
        "integrated_energy_j\n") < 0)
        return 0;
    AlcVec3 perpendicular = transverse_direction(alc_unit(c->direction));
    AlcVec3 wall = {perpendicular.x * c->radius_m,
                    perpendicular.y * c->radius_m,
                    perpendicular.z * c->radius_m};
    size_t count = alc_sample_count(c);
    for (size_t i = 0; i < count; ++i) {
        AlcState state = alc_state_at(c, alc_sample_time(c, i));
        AlcField field = alc_field_at(c, state.throttle, wall);
        if (fprintf(file, "%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,1,%.17g,%.17g\n",
            state.time_s, state.throttle, c->warp_factor, state.speed_m_s,
            state.speed_m_s / ALC_C, state.position_m.x, state.position_m.y,
            state.position_m.z, field.energy_density_j_m3,
            peak.energy_j * state.throttle * state.throttle) < 0)
            return 0;
    }
    return !ferror(file);
}

static int summary(FILE *file, const Options *o, AlcEnergy peak)
{
    const AlcConfig *c = &o->config;
    AlcState final = alc_state_at(c, c->duration_s);
    AlcVec3 n = alc_unit(c->direction);
    return fprintf(file,
        "{\n"
        "  \"schema_version\": 1,\n"
        "  \"simulator_version\": \"" ALC_VERSION "\",\n"
        "  \"model\": \"prescribed_alcubierre_metric\",\n"
        "  \"config\": {\n"
        "    \"radius_m\": %.17g, \"sigma_per_m\": %.17g, \"warp_factor\": %.17g,\n"
        "    \"direction\": [%.17g, %.17g, %.17g],\n"
        "    \"ramp_up_s\": %.17g, \"cruise_s\": %.17g, \"ramp_down_s\": %.17g,\n"
        "    \"duration_s\": %.17g, \"sample_dt_s\": %.17g, \"energy_intervals\": %u\n"
        "  },\n"
        "  \"samples\": %zu,\n"
        "  \"final_position_m\": [%.17g, %.17g, %.17g],\n"
        "  \"final_velocity_m_s\": %.17g,\n"
        "  \"full_throttle_energy_j\": %.17g,\n"
        "  \"full_throttle_mass_equivalent_kg\": %.17g,\n"
        "  \"quadrature\": {\n"
        "    \"method\": \"composite_simpson_wall_coordinate\",\n"
        "    \"fine_intervals\": %u, \"relative_difference\": %.17g,\n"
        "    \"converged\": %s, \"radial_min_m\": %.17g, \"radial_max_m\": %.17g\n"
        "  },\n"
        "  \"slice\": {\"enabled\": %s, \"time_s\": %.17g, \"points_per_axis\": %u,"
        " \"extent_m\": %.17g}\n"
        "}\n",
        c->radius_m, c->sigma_per_m, c->warp_factor,
        n.x, n.y, n.z, c->ramp_up_s, c->cruise_s, c->ramp_down_s,
        c->duration_s, c->sample_dt_s, c->energy_intervals, alc_sample_count(c),
        final.position_m.x, final.position_m.y, final.position_m.z,
        final.speed_m_s, peak.energy_j, peak.energy_j / (ALC_C * ALC_C),
        peak.intervals, peak.relative_difference, peak.converged ? "true" : "false",
        peak.radial_min_m, peak.radial_max_m, o->slice_path ? "true" : "false",
        o->slice_time_s, o->slice_points, o->slice_extent_m) >= 0 && !ferror(file);
}

static int slice(FILE *file, const Options *o)
{
    const AlcConfig *c = &o->config;
    AlcVec3 n = alc_unit(c->direction), p = transverse_direction(n);
    double throttle = alc_throttle(c, o->slice_time_s);
    if (fprintf(file, "longitudinal_m,transverse_m,shape,energy_density_j_m3,"
        "expansion_per_s,g_tt_m2_s2,shift_parallel_m_s\n") < 0)
        return 0;
    for (unsigned j = 0; j < o->slice_points; ++j) {
        double y = o->slice_extent_m * (2.0 * j / (o->slice_points - 1) - 1.0);
        for (unsigned i = 0; i < o->slice_points; ++i) {
            double x = o->slice_extent_m * (2.0 * i / (o->slice_points - 1) - 1.0);
            AlcVec3 offset = {x * n.x + y * p.x,
                              x * n.y + y * p.y,
                              x * n.z + y * p.z};
            AlcField field = alc_field_at(c, throttle, offset);
            if (fprintf(file, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                x, y, field.shape, field.energy_density_j_m3,
                field.expansion_per_s, field.metric[0][0],
                -ALC_C * c->warp_factor * throttle * field.shape) < 0)
                return 0;
        }
    }
    return !ferror(file);
}

int main(int argc, char **argv)
{
    Options o = {0};
    o.config = alc_default_config();
    o.output_path = "warp_telemetry.csv";
    o.summary_path = "warp_summary.json";
    o.slice_time_s = NAN;
    o.slice_extent_m = NAN;
    o.slice_points = 101;
    int parsed = parse(argc, argv, &o);
    if (parsed != 1)
        return parsed == 2 ? EXIT_SUCCESS : EXIT_FAILURE;
    char error[200];
    if (!alc_validate(&o.config, error, sizeof(error))) {
        fprintf(stderr, "configuration error: %s\n", error);
        return EXIT_FAILURE;
    }
    if (isnan(o.slice_time_s))
        o.slice_time_s = fmin(o.config.duration_s, o.config.ramp_up_s + 0.5 * o.config.cruise_s);
    if (isnan(o.slice_extent_m))
        o.slice_extent_m = o.config.radius_m + 6.0 / o.config.sigma_per_m;
    if ((o.slice_options && !o.slice_path) || o.slice_time_s < 0.0 ||
        o.slice_time_s > o.config.duration_s || o.slice_points < 3 ||
        o.slice_points > 501 || o.slice_extent_m <= 0.0 || o.slice_extent_m > 1e17) {
        fputs("invalid slice options: require --slice, time within run, 3..501 points,"
              " and extent in (0, 1e17] m\n", stderr);
        return EXIT_FAILURE;
    }
    const char *paths[3] = {o.output_path, o.summary_path, o.slice_path};
    size_t file_count = o.slice_path ? 3 : 2;
    for (size_t i = 0; i < file_count; ++i) {
        if (paths[i][0] == '\0' || strcmp(paths[i], "-") == 0) {
            fputs("outputs require nonempty file paths (stdout is not supported)\n", stderr);
            return EXIT_FAILURE;
        }
        for (size_t j = 0; j < i; ++j) {
            if (strcmp(paths[i], paths[j]) == 0) {
                fputs("output paths must be distinct\n", stderr);
                return EXIT_FAILURE;
            }
        }
    }
    AlcEnergy peak = alc_energy_at(&o.config, 1.0);
    if (!peak.converged) {
        fprintf(stderr, "energy quadrature did not converge (relative difference %.6g);"
                        " increase --energy-intervals\n", peak.relative_difference);
        return EXIT_FAILURE;
    }
    FILE *files[3] = {NULL, NULL, NULL};
    int created[3] = {0, 0, 0};
    int ok = 1;
    for (size_t i = 0; i < file_count; ++i) {
        files[i] = fopen(paths[i], "wx");
        if (files[i] == NULL) {
            perror(paths[i]);
            ok = 0;
            break;
        }
        created[i] = 1;
    }
    if (ok)
        ok = telemetry(files[0], &o.config, peak) && summary(files[1], &o, peak) &&
             (!o.slice_path || slice(files[2], &o));
    for (size_t i = 0; i < file_count; ++i) {
        if (files[i] != NULL && fclose(files[i]) != 0) {
            perror(paths[i]);
            ok = 0;
        }
    }
    if (!ok) {
        fputs("output failed; discarding files created by this run\n", stderr);
        for (size_t i = 0; i < file_count; ++i) {
            if (created[i] && remove(paths[i]) != 0)
                perror(paths[i]);
        }
        return EXIT_FAILURE;
    }
    AlcState final = alc_state_at(&o.config, o.config.duration_s);
    printf("UNSC Alcubierre Simulator %s | mathematical spacetime model\n"
           "Samples: %zu | final time: %.9g s\n"
           "Final position: [%.12g, %.12g, %.12g] m\n"
           "Full-throttle integrated Eulerian energy: %.12e J\n"
           "Quadrature relative difference: %.3e\n"
           "Telemetry: %s\nSummary: %s\n",
           ALC_VERSION, alc_sample_count(&o.config), final.time_s,
           final.position_m.x, final.position_m.y, final.position_m.z,
           peak.energy_j, peak.relative_difference, o.output_path, o.summary_path);
    if (o.slice_path)
        printf("Field slice: %s\n", o.slice_path);
    return EXIT_SUCCESS;
}
