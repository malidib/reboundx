/**
 * @file    stellar_evolution_sse.c
 * @brief   Very simplified stellar evolution update using analytic scalings.
 *
 * This operator updates a star's radius and luminosity based on its age
 * following approximate single-star evolution relations inspired by the
 * SSE formulation of \citet{Hurley2000}.  It is not a full
 * implementation of SSE, but provides time-dependent stellar properties
 * for testing wind mass loss.
 *
 * The main-sequence lifetime is approximated by
 *
 * .. math::
 *    t_{\rm MS} \approx 10^{10}\,\mathrm{yr}\,(M/M_\odot)^{-2.5}
 *
 * The zero-age main sequence radius and luminosity scale as
 *
 * .. math::
 *    R_{\rm ZAMS} \propto (M/M_\odot)^{0.8},\qquad
 *    L_{\rm ZAMS} \propto (M/M_\odot)^{3.5}.
 *
 * As the star ages on the main sequence, both its radius and luminosity
 * increase linearly with the fractional age $\tau=t/t_{\rm MS}$ until
 * $\tau=1$, after which they are held fixed.  The updated radius and
 * luminosity are stored both in the particle's ``r`` field and in the
 * parameters ``swml_R`` and ``swml_L`` so that the wind mass-loss
 * operator can access them.
 *
 * **Particle Parameters**
 *
 * ============================ =========== ====================================
 * Field (C type)               Required    Description
 * ============================ =========== ====================================
 * sse_age (double)             Yes         Current stellar age
 * ============================ =========== ====================================
 *
 * **Operator Parameters**
 *
 * ============================ =========== ====================================
 * Field (C type)               Required    Description
 * ============================ =========== ====================================
 * sse_Rsun (double)            No          Solar radius in code units (default 1)
 * sse_Lsun (double)            No          Solar luminosity in code units (default 1)
 * sse_R_coeff (double)        No          Multiplicative factor for radius scaling (default 1)
 * sse_R_exp (double)          No          Mass exponent for radius scaling (default 0.8)
 * sse_L_coeff (double)        No          Multiplicative factor for luminosity scaling (default 1)
 * sse_L_exp (double)          No          Mass exponent for luminosity scaling (default 3.5)
 * ============================ =========== ====================================
 *
 * The stellar mass is read directly from the particle's ``m`` field and is
 * assumed to be expressed in solar masses.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_stellar_evolution_sse(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    struct rebx_extras* const rebx = sim->extras;
    const int N_real = sim->N - sim->N_var;

    double Rsun = 1.;
    double Lsun = 1.;
    double R_coeff = 1.;
    double R_exp = 0.8;
    double L_coeff = 1.;
    double L_exp = 3.5;

    const double* Rsun_ptr = rebx_get_param(rebx, operator->ap, "sse_Rsun");
    const double* Lsun_ptr = rebx_get_param(rebx, operator->ap, "sse_Lsun");
    const double* R_coeff_ptr  = rebx_get_param(rebx, operator->ap, "sse_R_coeff");
    const double* R_exp_ptr    = rebx_get_param(rebx, operator->ap, "sse_R_exp");
    const double* L_coeff_ptr  = rebx_get_param(rebx, operator->ap, "sse_L_coeff");
    const double* L_exp_ptr    = rebx_get_param(rebx, operator->ap, "sse_L_exp");
    if (Rsun_ptr) Rsun = *Rsun_ptr;
    if (Lsun_ptr) Lsun = *Lsun_ptr;
    if (R_coeff_ptr) R_coeff = *R_coeff_ptr;
    if (R_exp_ptr)   R_exp   = *R_exp_ptr;
    if (L_coeff_ptr) L_coeff = *L_coeff_ptr;
    if (L_exp_ptr)   L_exp   = *L_exp_ptr;

    for (int i=0; i<N_real; i++){
        struct reb_particle* const p = &sim->particles[i];
        double* age_ptr = rebx_get_param(rebx, p->ap, "sse_age");
        if (age_ptr == NULL){
            continue;
        }

        *age_ptr += dt;
        const double mass_ratio = p->m; // assume simulation masses are in M_sun
        const double Rzams = R_coeff * Rsun * pow(mass_ratio, R_exp);
        const double Lzams = L_coeff * Lsun * pow(mass_ratio, L_exp);

        const double tms = 1e10 * pow(mass_ratio, -2.5); // yrs
        double tau = *age_ptr / tms;
        if (tau > 1.) tau = 1.;
        double R = Rzams * (1. + 5.*tau);
        double L = Lzams * (1. + 10.*tau);

        p->r = R;
        rebx_set_param_double(rebx, &p->ap, "swml_R", R);
        rebx_set_param_double(rebx, &p->ap, "swml_L", L);
    }
}
