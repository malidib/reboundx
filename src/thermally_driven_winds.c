/**
 * @file    thermally_driven_winds.c
 * @brief   Isotropic mass loss from thermally driven stellar winds.
 *
 * Default scaling (Parker‑like)
 *   dM/dt = - C_th * η *
 *           (R/R☉)^{α_R} (T/T☉)^{α_T} (M☉/M)^{α_M}  [M☉ yr⁻¹] .
 *
 * Operator‑level parameters
 * -------------------------
 *  tdw_const      (double)  – C_th, prefactor in M☉ yr⁻¹            (def: 2e‑14)
 *  tdw_Msun       (double)  – solar mass    in code units           (def: 1)
 *  tdw_Rsun       (double)  – solar radius  in code units           (def: 1)
 *  tdw_Tsun       (double)  – reference corona T in code units      (def: 1.5e6)
 *  tdw_year       (double)  – Julian year in code‑time units        (def: 1)
 *  tdw_alpha_R    (double)  – exponent α_R                          (def: 2)
 *  tdw_alpha_T    (double)  – exponent α_T                          (def: 1.5)
 *  tdw_alpha_M    (double)  – exponent α_M                          (def: 1)
 *  tdw_max_dlnM   (double)  – max |ΔM|/M per call (safety limiter)   (def: 0.1)
 *
 * Particle‑level parameters
 * -------------------------
 *  tdw_eta        (double, req.) – heating efficiency η
 *  tdw_T          (double, req.) – coronal temperature  (same units as Tsun)
 *  tdw_R          (double, req.) – stellar radius       (same units as Rsun)
 *
 * Notes
 * -----
 * * Mass is removed isotropically; no linear‑momentum recoil is applied.
 * * Operator is unit‑agnostic if you set the scaling constants consistently.
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"

/* Small helper: safe power with non‑positive base → returns 0 */
static inline double pow_safe(double x, double a){
    return (x > 0.0 && isfinite(x)) ? pow(x, a) : 0.0;
}

void rebx_thermally_driven_winds(struct reb_simulation* const sim,
                                 struct rebx_operator*  const op,
                                 const double                 dt)
{
    if (dt <= 0.0 || !isfinite(dt)) return;

    struct rebx_extras* const rx = sim->extras;
    const int N_real = sim->N - sim->N_var;                 /* ignore virtuals */

    /* ---------------- operator‑level parameters ------------------------ */
    double Msun      = 1.0;
    double Rsun      = 1.0;
    double Tsun      = 1.5e6;          /* K – typical solar coronal T      */
    double C0        = 2.0e-14;        /* M☉ yr⁻¹  → matches solar wind    */
    double year_len  = 1.0;

    double aR = 2.0;                   /* Parker: surface area             */
    double aT = 1.5;                   /* Parker: ρ_b c_s ∝ T^{3/2}        */
    double aM = 1.0;                   /* escape‑speed scaling             */
    double max_dlnM = 0.1;             /* safety: 10 % per call            */

    /* Retrieve overrides if present */
    #define GET(name,var) do{                              \
        const double* _p = rebx_get_param(rx, op->ap, name);\
        if(_p && isfinite(*_p) && *_p > 0.0) var = *_p;    \
    } while(0)

    GET("tdw_Msun"    , Msun     );
    GET("tdw_Rsun"    , Rsun     );
    GET("tdw_Tsun"    , Tsun     );
    GET("tdw_const"   , C0       );
    GET("tdw_year"    , year_len );
    GET("tdw_alpha_R" , aR       );
    GET("tdw_alpha_T" , aT       );
    GET("tdw_alpha_M" , aM       );
    GET("tdw_max_dlnM", max_dlnM );
    #undef GET

    /* Convert the prefactor to code‑mass per code‑time */
    const double pref = (C0 * Msun) / year_len;

    /* ---------------- iterate over real particles ---------------------- */
    for (int i = 0; i < N_real; i++){
        struct reb_particle* const p = &sim->particles[i];

        const double* eta_ptr = rebx_get_param(rx, p->ap, "tdw_eta");
        const double* T_ptr   = rebx_get_param(rx, p->ap, "tdw_T");
        const double* R_ptr   = rebx_get_param(rx, p->ap, "tdw_R");
        if (!eta_ptr || !T_ptr || !R_ptr) continue;

        const double eta = *eta_ptr;
        const double T   = *T_ptr;
        const double R   = *R_ptr;
        if (eta <= 0.0 || !isfinite(eta)) continue;
        if (T   <= 0.0 || !isfinite(T  )) continue;
        if (R   <= 0.0 || !isfinite(R  )) continue;
        if (p->m <= 0.0 || !isfinite(p->m)) continue;

        /* Dimensionless ratios and powers */
        const double RR = pow_safe(R / Rsun, aR);
        const double TT = pow_safe(T / Tsun, aT);
        const double MM = pow_safe(Msun / p->m, aM);

        /* Mass‑loss rate (code units) */
        const double mdot = -pref * eta * RR * TT * MM;

        /* Safety limiter */
        const double dM_lim = -max_dlnM * p->m;            /* negative value   */
        double dM = mdot * dt;
        if (dM < dM_lim) dM = dM_lim;

        /* Apply mass loss */
        p->m += dM;
        if (p->m < 0.0) p->m = 0.0;                       /* shouldn’t happen */
    }

    /* Recentre on centre of mass after mass has changed */
    reb_simulation_move_to_com(sim);
}
