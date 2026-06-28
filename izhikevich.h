#ifndef IZHIKEVICH_H
#define IZHIKEVICH_H

#include <stddef.h>
#include "plot_widget.h"

pointD *generate_sine_cos(size_t *out_n, pointD *out_v[3]);
pointD *generate_poly(size_t *out_n, pointD *out_v[3]);
pointD *generate_damped(size_t *out_n, pointD *out_v[3]);
pointD *generate_trig(size_t *out_n, pointD *out_v[3]);

#endif
