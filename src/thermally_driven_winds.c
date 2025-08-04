/**
 * @file    thermally_driven_winds.c
 * @brief   Isotropic mass loss from thermally driven stellar winds.
 *
 * Mass-loss rate:
 *   dM/dt = - C_th * eta * (R/R☉)^2 (T/T☉)^2 (M☉/M)  M☉ yr⁻¹
 *
 * Operator-level parameters
 * -------------------------
 *  tdw_const  (double, optional)  – prefactor in M☉ yr⁻¹ (default 1e-14)
 *  tdw_Msun   (double, optional)  – Solar mass  in code-mass units (default 1)
 *  tdw_Rsun   (double, optional)  – Solar radius in code-length units (default 1)
 *  tdw_Tsun   (double, optional)  – Solar temperature in code-temp units (default 1)
 *  tdw_year   (double, optional)  – Length of one Julian year in code-time
 *                                    units (default 1).
 *
 * Particle-level parameters
 * -------------------------
 *  tdw_eta    (double, required)  – dimensionless efficiency η
 *  tdw_T      (double, required)  – surface temperature (same units as Tsun)
 *  tdw_R      (double, required)  – stellar radius      (same units as Rsun)
 *
 * Notes
 * -----
 * * The operator is unit agnostic as long as the scaling parameters are set
 *   consistently.
 * * Mass is removed isotropically; no linear-momentum recoil is applied.
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"

/* ------------------------------------------------------------------------- */
/* Operator kernel                                                           */
/* ------------------------------------------------------------------------- */
void rebx_thermally_driven_winds(struct reb_simulation* const sim,
                                 struct rebx_operator*  const op,
                                 const double                 dt)
{
    if (dt <= 0.0 || !isfinite(dt)) return;

    struct rebx_extras* const rx = sim->extras;
    const int N_real = sim->N - sim->N_var;          /* ignore virtuals       */

    /* --------------------- operator-level parameters -------------------- */
    double Msun     = 1.0;    /* code-mass units per solar mass             */
    double Rsun     = 1.0;    /* code-length units per solar radius         */
    double Tsun     = 1.0;    /* code-temp   units per solar temperature    */
    double C0       = 1.0e-14;/* base prefactor [M☉ yr⁻¹]                    */
    double year_len = 1.0;    /* code-time units per Julian year            */

    const double* pM  = rebx_get_param(rx, op->ap, "tdw_Msun");
    const double* pR  = rebx_get_param(rx, op->ap, "tdw_Rsun");
    const double* pT  = rebx_get_param(rx, op->ap, "tdw_Tsun");
    const double* pC  = rebx_get_param(rx, op->ap, "tdw_const");
    const double* pYr = rebx_get_param(rx, op->ap, "tdw_year");

    if (pM  && isfinite(*pM ) && *pM  > 0.0) Msun     = *pM;
    if (pR  && isfinite(*pR ) && *pR  > 0.0) Rsun     = *pR;
    if (pT  && isfinite(*pT ) && *pT  > 0.0) Tsun     = *pT;
    if (pC  && isfinite(*pC ) && *pC  > 0.0) C0       = *pC;
    if (pYr && isfinite(*pYr) && *pYr > 0.0) year_len = *pYr;

    /* Prefactor converted to code-mass per code-time */
    const double pref = (C0 * Msun) / year_len;      /* [code-mass / code-time] */

    /* --------------------- loop over genuine particles ------------------ */
    for (int i = 0; i < N_real; i++){
        struct reb_particle* const p = &sim->particles[i];

        /* Fetch particle parameters */
        const double* eta_ptr = rebx_get_param(rx, p->ap, "tdw_eta");
        const double* T_ptr   = rebx_get_param(rx, p->ap, "tdw_T");
        const double* R_ptr   = rebx_get_param(rx, p->ap, "tdw_R");

        if (!eta_ptr || !T_ptr || !R_ptr) continue;

        const double eta = *eta_ptr;
        const double T   = *T_ptr;
        const double R   = *R_ptr;

        /* Sanity checks */
        if (!isfinite(eta) || eta <= 0.0) continue;
        if (!isfinite(T)   || T   <= 0.0) continue;
        if (!isfinite(R)   || R   <= 0.0) continue;
        if (p->m <= 0.0 || !isfinite(p->m)) continue;

        /* Thermal wind mass-loss rate in code units */
        const double mdot = -pref * eta
                            * (R / Rsun) * (R / Rsun)
                            * (T / Tsun) * (T / Tsun)
                            * (Msun / p->m);     /* [code-mass / code-time] */

        const double dM = mdot * dt;             /* negative number         */

        /* Prevent negative masses */
        if (p->m + dM <= 0.0) {
            p->m = 0.0;
        } else {
            p->m += dM;
        }
    }

    /* Recentre: mass distribution changed */
    reb_simulation_move_to_com(sim);
}

