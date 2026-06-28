#include "izhikevich.h"
#include <stdlib.h>
#include <math.h>

#define NPOINTS 256

pointD *generate_sine_cos(size_t *out_n, pointD *out_v[3])
{
    int n = NPOINTS;
    out_v[0] = malloc(sizeof(pointD) * n);
    out_v[1] = malloc(sizeof(pointD) * n);
    out_v[2] = malloc(sizeof(pointD) * n);
    for (int i = 0; i < n; ++i) {
        double x = 4.0 * M_PI * i / (n - 1);
        out_v[0][i] = (pointD){ x, sin(x) };
        out_v[1][i] = (pointD){ x, cos(x) };
        out_v[2][i] = (pointD){ x, sin(x) * cos(x) };
    }
    *out_n = n;
    return out_v[1];
}

pointD *generate_poly(size_t *out_n, pointD *out_v[3])
{
    int n = NPOINTS;
    out_v[0] = malloc(sizeof(pointD) * n);
    out_v[1] = malloc(sizeof(pointD) * n);
    out_v[2] = malloc(sizeof(pointD) * n);
    for (int i = 0; i < n; ++i) {
        double x = -5.0 + 10.0 * i / (n - 1);
        out_v[0][i] = (pointD){ x, x * x - 4.0 };
        out_v[1][i] = (pointD){ x, x * x * x / 10.0 };
        out_v[2][i] = (pointD){ x, x };
    }
    *out_n = n;
    return out_v[1];
}

pointD *generate_damped(size_t *out_n, pointD *out_v[3])
{
    int n = NPOINTS;
    out_v[0] = malloc(sizeof(pointD) * n);
    out_v[1] = malloc(sizeof(pointD) * n);
    out_v[2] = malloc(sizeof(pointD) * n);
    for (int i = 0; i < n; ++i) {
        double x = 4.0 * M_PI * i / (n - 1);
        out_v[0][i] = (pointD){ x, exp(-x / 5.0) * sin(x) };
        out_v[1][i] = (pointD){ x, exp(-x / 10.0) * cos(x) };
        out_v[2][i] = (pointD){ x, exp(-x / 7.0) * sin(2.0 * x) };
    }
    *out_n = n;
    return out_v[1];
}

pointD *generate_trig(size_t *out_n, pointD *out_v[3])
{
    int n = NPOINTS;
    out_v[0] = malloc(sizeof(pointD) * n);
    out_v[1] = malloc(sizeof(pointD) * n);
    out_v[2] = malloc(sizeof(pointD) * n);
    for (int i = 0; i < n; ++i) {
        double x = 4.0 * M_PI * i / (n - 1);
        out_v[0][i] = (pointD){ x, sin(2.0 * x) };
        out_v[1][i] = (pointD){ x, sin(3.0 * x) };
        out_v[2][i] = (pointD){ x, cos(x) * sin(2.0 * x) };
    }
    *out_n = n;
    return out_v[1];
}
