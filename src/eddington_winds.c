/**
 * @file    eddington_winds.c
 * @brief   Stellar mass loss driven by super-Eddington luminosities.
 *
 * Mass-loss rate:
 *   dM/dt = - C_Edd * max(0, L/L_Edd - 1)   [M☉ yr⁻¹].
 *
 * L_Edd is approximated as
 *   L_Edd = coeff * (M/M☉) L☉.
 *
 * Operator-level parameters
 * -------------------------
 *  eddw_const       (double) – C_Edd prefactor in M☉ yr⁻¹        (def: 1e-5)
 *  eddw_Msun        (double) – solar mass in code units          (def: 1)
 *  eddw_Lsun        (double) – solar luminosity in code units    (def: 1)
 *  eddw_LEdd_coeff  (double) – coeff for L_Edd in L☉            (def: 3.2e4)
 *  eddw_year        (double) – Julian year in code-time units    (def: 1)
 *  eddw_max_dlnM    (double) – max |ΔM|/M per call               (def: 0.1)
 *
 * Particle-level parameters
 * ------------------------
 *  eddw_L           (double, req.) – stellar luminosity (same units as Lsun)
 *
 * Notes
 * -----
 * * Mass is removed isotropically; no linear-momentum recoil is applied.
 * * Operator is unit-agnostic if scaling constants are set consistently.
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_eddington_winds(struct reb_simulation* const sim,
                          struct rebx_operator*  const op,
                          const double                 dt)
{
    if (dt <= 0.0 || !isfinite(dt)) return;

    struct rebx_extras* const rx = sim->extras;
    const int N_real = sim->N - sim->N_var; /* ignore virtuals */

    /* Operator-level defaults */
    double Msun       = 1.0;
    double Lsun       = 1.0;
    double C0         = 1e-5;      /* M☉ yr⁻¹ */
    double LEdd_coeff = 3.2e4;     /* in L☉ */
    double year_len   = 1.0;
    double max_dlnM   = 0.1;

    /* Retrieve overrides if present */
    #define GET(name,var) do{                              \
        const double* _p = rebx_get_param(rx, op->ap, name);\
        if(_p && isfinite(*_p) && *_p > 0.0) var = *_p;    \
    } while(0)

    GET("eddw_Msun"      , Msun      );
    GET("eddw_Lsun"      , Lsun      );
    GET("eddw_const"     , C0        );
    GET("eddw_LEdd_coeff", LEdd_coeff);    /* dimensionless */
    GET("eddw_year"      , year_len  );
    GET("eddw_max_dlnM"  , max_dlnM  );
    #undef GET

    const double pref = (C0 * Msun) / year_len; /* convert to code-mass/time */

    for (int i = 0; i < N_real; i++){
        struct reb_particle* const p = &sim->particles[i];

        const double* L_ptr = rebx_get_param(rx, p->ap, "eddw_L");
        if (!L_ptr) continue;
        const double L = *L_ptr;
        if (!isfinite(L) || L <= 0.0) continue;
        if (p->m <= 0.0 || !isfinite(p->m)) continue;

        const double L_Edd = LEdd_coeff * (p->m / Msun) * Lsun;
        if (!isfinite(L_Edd) || L_Edd <= 0.0) continue;
        if (L <= L_Edd) continue; /* below Eddington, no wind */

        const double Gamma = L / L_Edd;
        const double mdot  = -pref * (Gamma - 1.0);

        double dM = mdot * dt;                 /* negative */
        const double dM_lim = -max_dlnM * p->m;
        if (dM < dM_lim) dM = dM_lim;

        p->m += dM;
        if (p->m < 0.0) p->m = 0.0;
    }

    /* Recentre after mass change */
    reb_simulation_move_to_com(sim);
}

