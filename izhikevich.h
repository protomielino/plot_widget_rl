#ifndef IZHIKEVICH_H
#define IZHIKEVICH_H

#include <stddef.h>
#include "plot_widget.h"

typedef enum {
    IZK_MODEL_TONIC_SPIKING,
    IZK_MODEL_PHASIC_SPIKING,
    IZK_MODEL_TONIC_BURSTING,
    IZK_MODEL_PHASIC_BURSTING,
    IZK_MODEL_MIXED_MODE,
    IZK_MODEL_SPIKE_FREQUENCY_ADAPTATION,
    IZK_MODEL_CLASS_1,
    IZK_MODEL_CLASS_2,
    IZK_MODEL_SPIKE_LATENCY,
    IZK_MODEL_SUBTHRESHOLD_OSCILLATIONS,
    IZK_MODEL_RESONATOR,
    IZK_MODEL_INTEGRATOR,
    IZK_MODEL_REBOUND_SPIKE,
    IZK_MODEL_REBOUND_BURST,
    IZK_MODEL_THRESHOLD_VARIABILITY,
    IZK_MODEL_BISTABILITY,
    IZK_MODEL_DAP,
    IZK_MODEL_ACCOMODATION,
    IZK_MODEL_INHIBITION_INDUCED_SPIKING,
    IZK_MODEL_INHIBITION_INDUCED_BURSTING,
    IZK_MODEL_N
} izk_model_e;

/* Each generator fills out_v[3] with time-series for (I, V, u)
 * and returns out_v[1] (V) as a convenience.
 * *out_n is set to the number of points.
 */
pointD *generate_izk_tonic_spike(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_phasic_spike(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_tonic_bursting(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_phasic_bursting(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_mixed_mode(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_spike_freq_adapt(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_class1_exc(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_class2_exc(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_spike_latency(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_subthr_osc(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_resonator(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_integrator(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_rebound_spike(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_rebound_burst(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_thresh_variability(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_bistability(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_DAP(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_accomodation(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_inh_induced_sp(size_t *out_n, pointD *out_v[3]);
pointD *generate_izk_inh_induced_brst(size_t *out_n, pointD *out_v[3]);

#endif /* IZHIKEVICH_H */
