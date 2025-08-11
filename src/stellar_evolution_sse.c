/**
 * @file    stellar_evolution_sse.c
 * @brief   Very simplified stellar evolution update using analytic scalings.
 *
 * This operator updates a star's radius and luminosity using simple
 * power-law scalings with its mass, inspired by the SSE formulation of
 * \citet{Hurley2000}.  It is not a full implementation of SSE, but
 * provides mass-dependent stellar properties for testing wind mass loss.
 *
 * The zero-age main sequence radius and luminosity scale approximately as
 *
 * 
 *    R \propto (M/M_\odot)^{0.8},\qquad
 *    L \propto (M/M_\odot)^{3.5}.
 *
 * ============================ =========== ====================================
 * Field (C type)               Required    Description
 * ============================ =========== ====================================
 * sse_Msun (double)           No          Solar mass in code units (default 1)
 * sse_Rsun (double)           No          Solar radius in code units (default 1)
 * sse_Lsun (double)           No          Solar luminosity in code units (default 1)
 * sse_R_coeff (double)        No          Multiplicative factor for radius scaling (default 1)
 * sse_R_exp (double)          No          Mass exponent for radius scaling (default 0.8)
 * sse_L_coeff (double)        No          Multiplicative factor for luminosity scaling (default 1)
 * sse_L_exp (double)          No          Mass exponent for luminosity scaling (default 3.5)
 * ============================ =========== ====================================
 *
 * The stellar mass is read directly from the particle's ``m`` field and is
 * converted from code units using ``sse_Msun``; users therefore do not need
 * to work in solar‑mass units.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_stellar_evolution_sse(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    struct rebx_extras* const rebx = sim->extras;
    const int N_real = sim->N - sim->N_var;

    double Msun = 1.;
    double Rsun = 1.;
    double Lsun = 1.;
    // Default scaling coefficients/exponents
    const double R_coeff_default = 1.;
    const double R_exp_default   = 0.8;
    const double L_coeff_default = 1.;
    const double L_exp_default   = 3.5;

    const double* Msun_ptr = rebx_get_param(rebx, operator->ap, "sse_Msun");
    const double* Rsun_ptr = rebx_get_param(rebx, operator->ap, "sse_Rsun");
    const double* Lsun_ptr = rebx_get_param(rebx, operator->ap, "sse_Lsun");
    if (Msun_ptr) Msun = *Msun_ptr;
    if (Rsun_ptr) Rsun = *Rsun_ptr;
    if (Lsun_ptr) Lsun = *Lsun_ptr;

    for (int i=0; i<N_real; i++){
        struct reb_particle* const p = &sim->particles[i];
        double R_coeff = R_coeff_default;
        double R_exp   = R_exp_default;
        double L_coeff = L_coeff_default;
        double L_exp   = L_exp_default;
        const double* R_coeff_ptr = rebx_get_param(rebx, p->ap, "sse_R_coeff");
        const double* R_exp_ptr   = rebx_get_param(rebx, p->ap, "sse_R_exp");
        const double* L_coeff_ptr = rebx_get_param(rebx, p->ap, "sse_L_coeff");
        const double* L_exp_ptr   = rebx_get_param(rebx, p->ap, "sse_L_exp");
        if (R_coeff_ptr) R_coeff = *R_coeff_ptr;
        if (R_exp_ptr)   R_exp   = *R_exp_ptr;
        if (L_coeff_ptr) L_coeff = *L_coeff_ptr;
        if (L_exp_ptr)   L_exp   = *L_exp_ptr;
        const double mass_ratio = p->m / Msun;
        double R = R_coeff * Rsun * pow(mass_ratio, R_exp);
        double L = L_coeff * Lsun * pow(mass_ratio, L_exp);
        p->r = R;
        rebx_set_param_double(rebx, &p->ap, "sse_L", L);
    }
}
