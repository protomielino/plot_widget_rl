#include "izhikevich.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ================================================================== */
/*  (A) tonic spiking                                                  */
/* ================================================================== */
pointD *generate_izk_tonic_spike(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.02, b = 0.2, c = -65.0, d = 6.0;
    double V = -70.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 100.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); /* I */
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); /* V */
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); /* u */
    const double T1 = t_end / 10.0;

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double I = (i * tau > T1) ? 14.0 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (B) phasic spiking                                                 */
/* ================================================================== */
pointD *generate_izk_phasic_spike(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.02, b = 0.25, c = -65.0, d = 6.0;
    double V = -64.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 200.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = 20.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double I = (i * tau > T1) ? 0.5 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (C) tonic bursting                                                 */
/* ================================================================== */
pointD *generate_izk_tonic_bursting(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.02, b = 0.2, c = -50.0, d = 2.0;
    double V = -70.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 220.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = 22.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double I = (i * tau > T1) ? 15.0 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (D) phasic bursting                                                */
/* ================================================================== */
pointD *generate_izk_phasic_bursting(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.02, b = 0.25, c = -55.0, d = 0.05;
    double V = -64.0, u = b * V;
    const double tau = 0.2, t_start = 0.0, t_end = 200.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = 20.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double I = (i * tau > T1) ? 0.6 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (E) mixed mode                                                     */
/* ================================================================== */
pointD *generate_izk_mixed_mode(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.02, b = 0.2, c = -55.0, d = 4.0;
    double V = -70.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 160.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = t_end / 10.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double I = (i * tau > T1) ? 10.0 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (F) spike frequency adaptation                                     */
/* ================================================================== */
pointD *generate_izk_spike_freq_adapt(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.01, b = 0.2, c = -65.0, d = 8.0;
    double V = -70.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 85.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = t_end / 10.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double I = (i * tau > T1) ? 30.0 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (G) Class 1 excitable                                              */
/* ================================================================== */
pointD *generate_izk_class1_exc(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.02, b = -0.1, c = -55.0, d = 6.0;
    double V = -60.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 300.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = 30.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = (t > T1) ? 0.075 * (t - T1) : 0.0;
        V = V + tau * (0.04 * V * V + 4.1 * V + 108.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (H) Class 2 excitable                                              */
/* ================================================================== */
pointD *generate_izk_class2_exc(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.2, b = 0.26, c = -65.0, d = 0.0;
    double V = -64.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 300.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = 30.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = (t > T1) ? -0.5 + 0.015 * (t - T1) : -0.5;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (I) spike latency                                                  */
/* ================================================================== */
pointD *generate_izk_spike_latency(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.02, b = 0.2, c = -65.0, d = 6.0;
    double V = -70.0, u = b * V;
    const double tau = 0.2, t_start = 0.0, t_end = 100.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = t_end / 10.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = ((t > T1) && (t < T1 + 3.0)) ? 7.04 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (J) subthreshold oscillations                                      */
/* ================================================================== */
pointD *generate_izk_subthr_osc(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.05, b = 0.26, c = -60.0, d = 0.0;
    double V = -62.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 200.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = t_end / 10.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = ((t > T1) && (t < T1 + 5.0)) ? 2.0 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (K) resonator                                                      */
/* ================================================================== */
pointD *generate_izk_resonator(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.1, b = 0.26, c = -60.0, d = -1.0;
    double V = -62.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 400.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = t_end / 10.0;
    const double T2 = T1 + 20.0;
    const double T3 = 0.7 * t_end;
    const double T4 = T3 + 40.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = (((t > T1) && (t < T1 + 4.0)) ||
                    ((t > T2) && (t < T2 + 4.0)) ||
                    ((t > T3) && (t < T3 + 4.0)) ||
                    ((t > T4) && (t < T4 + 4.0))) ? 0.65 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (L) integrator                                                     */
/* ================================================================== */
pointD *generate_izk_integrator(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.02, b = -0.1, c = -55.0, d = 6.0;
    double V = -60.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 100.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = t_end / 11.0;
    const double T2 = T1 + 5.0;
    const double T3 = 0.7 * t_end;
    const double T4 = T3 + 10.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = (((t > T1) && (t < T1 + 2.0)) ||
                    ((t > T2) && (t < T2 + 2.0)) ||
                    ((t > T3) && (t < T3 + 2.0)) ||
                    ((t > T4) && (t < T4 + 2.0))) ? 9.0 : 0.0;
        V = V + tau * (0.04 * V * V + 4.1 * V + 108.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (M) rebound spike                                                  */
/* ================================================================== */
pointD *generate_izk_rebound_spike(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.03, b = 0.25, c = -60.0, d = 4.0;
    double V = -64.0, u = b * V;
    const double tau = 0.2, t_start = 0.0, t_end = 200.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = 20.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = ((t > T1) && (t < T1 + 5.0)) ? -15.0 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (N) rebound burst                                                  */
/* ================================================================== */
pointD *generate_izk_rebound_burst(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.03, b = 0.25, c = -52.0, d = 0.0;
    double V = -64.0, u = b * V;
    const double tau = 0.2, t_start = 0.0, t_end = 200.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = 20.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = ((t > T1) && (t < T1 + 5.0)) ? -15.0 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (O) threshold variability                                          */
/* ================================================================== */
pointD *generate_izk_thresh_variability(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.03, b = 0.25, c = -60.0, d = 4.0;
    double V = -64.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 100.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = (((t > 10.0) && (t < 15.0)) || ((t > 80.0) && (t < 85.0))) ? 1.0 :
                   ((t > 70.0) && (t < 75.0)) ? -6.0 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (P) bistability                                                    */
/* ================================================================== */
pointD *generate_izk_bistability(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.1, b = 0.26, c = -60.0, d = 0.0;
    double V = -61.0, u = b * V;
    const double tau = 0.25, t_start = 0.0, t_end = 300.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = t_end / 8.0;
    const double T2 = 216.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = (((t > T1) && (t < T1 + 5.0)) || ((t > T2) && (t < T2 + 5.0))) ? 1.24 : 0.24;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (Q) DAP (depolarizing after-potential)                             */
/* ================================================================== */
pointD *generate_izk_DAP(size_t *out_n, pointD *out_v[3])
{
    const double a = 1.0, b = 0.2, c = -60.0, d = -21.0;
    double V = -70.0, u = b * V;
    const double tau = 0.1, t_start = 0.0, t_end = 50.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;
    const double T1 = 10.0;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = (fabs(t - T1) < 1.0) ? 20.0 : 0.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (R) accommodation                                                  */
/* ================================================================== */
pointD *generate_izk_accomodation(size_t *out_n, pointD *out_v[3])
{
    const double a = 0.02, b = 1.0, c = -55.0, d = 4.0;
    double V = -65.0, u = -16.0; /* note: u = -16, not b*V */
    const double tau = 0.5, t_start = 0.0, t_end = 400.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    double *II = malloc(sizeof(double) * steps);
    if (!VV || !uu || !II) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I;
        if (t < 200.0)         I = t / 25.0;
        else if (t < 300.0)    I = 0.0;
        else if (t < 312.5)    I = (t - 300.0) / 12.5 * 4.0;
        else                   I = 0.0;

        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * (V + 65.0));
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        II[i] = I;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu); free(II);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (S) inhibition-induced spiking                                     */
/* ================================================================== */
pointD *generate_izk_inh_induced_sp(size_t *out_n, pointD *out_v[3])
{
    const double a = -0.02, b = -1.0, c = -60.0, d = 8.0;
    double V = -63.8, u = b * V;
    const double tau = 0.5, t_start = 0.0, t_end = 350.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = ((t < 50.0) || (t > 250.0)) ? 80.0 : 75.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}

/* ================================================================== */
/*  (T) inhibition-induced bursting                                    */
/* ================================================================== */
pointD *generate_izk_inh_induced_brst(size_t *out_n, pointD *out_v[3])
{
    const double a = -0.026, b = -1.0, c = -45.0, d = -2.0;
    double V = -63.8, u = b * V;
    const double tau = 0.5, t_start = 0.0, t_end = 350.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps);
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps);

    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!VV || !uu) { fprintf(stderr, "Insufficient memory\n"); return NULL; }

    for (int i = 0; i < steps; ++i) {
        double t = i * tau;
        double I = ((t < 50.0) || (t > 250.0)) ? 80.0 : 75.0;
        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);
        if (V > 30.0) { VV[i] = 30.0; V = c; u = u + d; }
        else { VV[i] = V; }
        uu[i] = u;
        out_v[0][i] = (pointD){ (double)i, I };
        out_v[1][i] = (pointD){ (double)i, VV[i] };
        out_v[2][i] = (pointD){ (double)i, uu[i] };
    }
    free(VV); free(uu);
    *out_n = steps;
    return out_v[1];
}
