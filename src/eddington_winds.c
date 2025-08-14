/**
 * @file    eddington_winds.c
 * @brief   Stellar mass loss driven by super-Eddington luminosities.
 *
 * Mass-loss rate:
 *   dM/dt = - C_edd * max(0, L/L_Edd - 1)   [M☉ yr⁻¹]
 *
 * L_Edd ≈ coeff * (M/M☉) L☉  (electron-scattering Eddington limit).
 *
 * Operator-level parameters (canonical; synonyms accepted)
 * -------------------------------------------------------
 *  edw_const        (double) – C_edd in M☉/yr              (def: 1e-6)
 *      synonyms: eddw_const, ew_const
 *  edw_Msun         (double) – solar mass in code units    (def: 1)
 *      synonyms: eddw_Msun, ew_Msun
 *  edw_Lsun         (double) – solar luminosity in code    (def: 1)
 *      synonyms: eddw_Lsun, ew_Lsun
 *  edw_year         (double) – year length in code time    (def: 1)
 *      synonyms: eddw_year, ew_year
 *  edw_Ledd_coeff   (double) – L_Edd per Msun in L☉/M☉     (def: 3.2e4)
 *      synonyms: edw_Ledd0, edw_Ledd, eddw_LEdd_coeff, ew_Ledd_coeff
 *  edw_max_dlnM     (double) – max |ΔM|/M per call         (def: 0.1)
 *
 * Particle-level parameters
 * -------------------------
 *  sse_L           (double, req.) – stellar luminosity (units of L☉)
 *
 * Notes
 * -----
 * * Mass is removed isotropically; no linear-momentum recoil.
 * * Unit-agnostic if scaling constants are set consistently.
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

    /* -------- operator-level parameters (with synonyms) -------- */
    double Msun       = 1.0;
    double Lsun       = 1.0;
    double year_len   = 1.0;
    double C0         = 1e-6;     /* M☉/yr (conservative default) */
    double Ledd_coeff = 3.2e4;    /* L☉/M☉ */
    double max_dlnM   = 0.1;

    #define GET(name,var) do{ \
        const double* _p = rebx_get_param(rx, op->ap, name); \
        if (_p && isfinite(*_p) && *_p > 0.0) var = *_p; \
    } while(0)

    /* canonical + synonyms */
    GET("edw_Msun",       Msun      ); GET("eddw_Msun",       Msun      ); GET("ew_Msun",       Msun      );
    GET("edw_Lsun",       Lsun      ); GET("eddw_Lsun",       Lsun      ); GET("ew_Lsun",       Lsun      );
    GET("edw_year",       year_len  ); GET("eddw_year",       year_len  ); GET("ew_year",       year_len  );
    GET("edw_const",      C0        ); GET("eddw_const",      C0        ); GET("ew_const",      C0        );
    GET("edw_Ledd_coeff", Ledd_coeff); GET("edw_Ledd0",       Ledd_coeff); GET("edw_Ledd",      Ledd_coeff);
    GET("eddw_LEdd_coeff",Ledd_coeff); GET("ew_Ledd_coeff",   Ledd_coeff);
    GET("edw_max_dlnM",   max_dlnM  );

    #undef GET

    const double pref = (C0 * Msun) / year_len; /* code mass per code time */

    /* ----------------- iterate over real particles ----------------- */
    for (int i = 0; i < N_real; i++){
        struct reb_particle* const p = &sim->particles[i];
        if (p->m <= 0.0 || !isfinite(p->m)) continue;

        /* Stellar luminosity supplied via sse_L (e.g. from SSE operator) */
        const double* Lp = rebx_get_param(rx, p->ap, "sse_L");
        if (!Lp) continue;

        const double L = *Lp;
        if (!isfinite(L) || L <= 0.0) continue;

        const double Ledd = Ledd_coeff * (p->m / Msun) * Lsun;
        if (!isfinite(Ledd) || Ledd <= 0.0) continue;

        const double Gamma = L / Ledd;
        if (!(Gamma > 1.0)) continue; /* at/below Eddington => no wind */

        const double mdot = -pref * (Gamma - 1.0); /* Mdot ≤ 0 */
        double dM = mdot * dt;

        /* safety limiter */
        const double dM_lim = -max_dlnM * p->m;
        if (dM < dM_lim) dM = dM_lim;

        p->m += dM;
        if (p->m < 0.0) p->m = 0.0;
    }

    reb_simulation_move_to_com(sim);
}
