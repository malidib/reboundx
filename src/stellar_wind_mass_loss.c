/**
 * @file    stellar_wind_mass_loss.c
 * @brief   Continuous stellar wind mass loss using Reimers' law.
 *
 * This operator decreases a star's mass every timestep according
 * to the Reimers wind prescription
 *
 *     \dot M = -4\times10^{-13}\,\eta\,\frac{L}{L_\odot}\frac{R}{R_\odot}\frac{M_\odot}{M}\;M_\odot\,\mathrm{yr}^{-1}
 *
 * where the user supplies the dimensionless efficiency \eta as well
 * as the stellar luminosity and radius.  All quantities should be
 * given in the same units as used in the simulation.  The constant
 * factor may be overridden via the operator parameter ``swml_const``.
 *
 * **Particle Parameters**
 *
 * ============================ =========== ====================================
 * Field (C type)               Required    Description
 * ============================ =========== ====================================
 * swml_eta (double)           Yes         Wind efficiency parameter \eta
 * swml_L   (double)           Yes         Stellar luminosity
 * swml_R   (double)           Yes         Stellar radius
 * ============================ =========== ====================================
 *
 * **Operator Parameters**
 *
 * ============================ =========== ====================================
 * Field (C type)               Required    Description
 * ============================ =========== ====================================
 * swml_const (double)          No          Prefactor in Msun/yr (default 4e-13)
 * swml_Msun  (double)          No          Solar mass in code units (default 1)
 * swml_Rsun  (double)          No          Solar radius in code units (default 1)
 * swml_Lsun  (double)          No          Solar luminosity in code units (default 1)
 * ============================ =========== ====================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_stellar_wind_mass_loss(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    struct rebx_extras* const rebx = sim->extras;
    const int N_real = sim->N - sim->N_var;

    double Msun = 1.;
    double Rsun = 1.;
    double Lsun = 1.;
    double constant = 4e-13; // default Msun/yr

    const double* Msun_ptr = rebx_get_param(rebx, operator->ap, "swml_Msun");
    const double* Rsun_ptr = rebx_get_param(rebx, operator->ap, "swml_Rsun");
    const double* Lsun_ptr = rebx_get_param(rebx, operator->ap, "swml_Lsun");
    const double* const_ptr = rebx_get_param(rebx, operator->ap, "swml_const");
    if (Msun_ptr) Msun = *Msun_ptr;
    if (Rsun_ptr) Rsun = *Rsun_ptr;
    if (Lsun_ptr) Lsun = *Lsun_ptr;
    if (const_ptr) constant = *const_ptr;

    for (int i=0; i<N_real; i++){
        struct reb_particle* const p = &sim->particles[i];
        const double* eta_ptr = rebx_get_param(rebx, p->ap, "swml_eta");
        const double* L_ptr   = rebx_get_param(rebx, p->ap, "swml_L");
        const double* R_ptr   = rebx_get_param(rebx, p->ap, "swml_R");
        if (eta_ptr && L_ptr && R_ptr){
            if (p->m <= 0.) continue;
            const double eta = *eta_ptr;
            const double L   = *L_ptr;
            const double R   = *R_ptr;
            double mdot = -constant * eta * (L/Lsun) * (R/Rsun) / (p->m/Msun); // Msun/yr
            double dM = mdot * dt;
            if (p->m + dM <= 0.){
                dM = -p->m;
            }
            p->m += dM;
        }
    }
    reb_simulation_move_to_com(sim);
}
