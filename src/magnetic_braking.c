/**
 * @file    magnetic_braking.c
 * @brief   Verbunt–Zwaan / Kawaler magnetic‑braking torque with saturation.
 *
 * dJ/dt = -K R^{1/2} M^{-1/2} ×
 *         { Ω^3                     , Ω ≤ Ω_sat
 *         { Ω Ω_sat^2              , Ω > Ω_sat
 *
 * where K, R, M are in consistent units.
 *
 * Particle‑level parameters  (all scalars are stored as double in REBOUNDx)
 * ------------------------------------------------------------------------
 * mb_on          (double, required ≠0)  – enable braking on this particle
 * mb_convective  (double, required ≠0)  – star has a convective envelope
 * mb_omega_sat   (double, optional)     – saturation angular velocity
 * mb_tau_conv    (double, optional)     – convective turnover time
 * I              (double, required)     – moment of inertia
 * Omega          (reb_vec3d, required)  – spin angular‑velocity vector
 *
 * Operator‑level parameters
 * -------------------------
 * mb_K           (double, optional)     – braking constant in cgs units
 *                                             (default 2.7e47).
 * mb_Msun        (double, optional)     – solar mass  in code‑mass units
 *                                             (default 1).
 * mb_Rsun        (double, optional)     – solar radius in code‑length units
 *                                             (default 1).
 * mb_year        (double, optional)     – Julian year in code‑time units
 *                                             (default 1).
 * mb_Rossby_sat  (double, optional)     – critical Rossby number for
 *                                         saturation (default 0.1)
 *
 * All quantities are automatically converted to the simulation's unit
 * system; users need not manually rescale ``mb_K``.
*/

#include <math.h>
#include <stdio.h>
#include "rebound.h"
#include "reboundx.h"

static inline void apply_magnetic_brake(struct reb_particle*    const p,
                                        struct rebx_extras*     const rx,
                                        const double            K,
                                        const double            dt,
                                        const double            Rossby_sat)
{
    /* ----------------------- fetch flags & properties ------------------- */
    const double* mb_on  = rebx_get_param(rx, p->ap, "mb_on");
    if (!mb_on || *mb_on == 0.0) return;

    const double* conv   = rebx_get_param(rx, p->ap, "mb_convective");
    if (!conv || *conv == 0.0) return;

    double* I_ptr = rebx_get_param(rx, p->ap, "I");
    if (!I_ptr || *I_ptr <= 0.0 || !isfinite(*I_ptr)) return;

    struct reb_vec3d* Omega = rebx_get_param_vec(rx, p->ap, "Omega");
    if (!Omega) return;

    const double R = p->r;
    const double M = p->m;
    if (R <= 0.0 || M <= 0.0 || !isfinite(R) || !isfinite(M)) return;

    const double omega = sqrt(Omega->x*Omega->x + Omega->y*Omega->y + Omega->z*Omega->z);
    if (!isfinite(omega) || omega == 0.0) return;

    const double* sat_ptr   = rebx_get_param(rx, p->ap, "mb_omega_sat");
    const double* tau_ptr   = rebx_get_param(rx, p->ap, "mb_tau_conv");
    double omega_sat = INFINITY;
    if (tau_ptr && *tau_ptr > 0.0 && isfinite(*tau_ptr)){
        omega_sat = 2.0*M_PI/(Rossby_sat * (*tau_ptr));
    }
    if (sat_ptr && *sat_ptr > 0.0 && isfinite(*sat_ptr))
        omega_sat = *sat_ptr;

    /* --------------------------- torque magnitude ----------------------- */
    double omega_term = omega*omega*omega;          /* unsaturated default  */
    if (omega > omega_sat)                          /* saturated regime     */
        omega_term = omega * omega_sat * omega_sat;

    const double torque = -K * sqrt(R) * (1.0 / sqrt(M)) * omega_term;

    /* --------------------------- apply dΩ/dt ---------------------------- */
    const double domega_dt = torque / (*I_ptr);     /* scalar rate (negative) */
    const double scale     = 1.0 + (domega_dt / omega)*dt;

    /* Prevent overshoot that could flip the spin */
    if (scale <= 0.0){
        fprintf(stderr, "[magnetic_braking] spin update skipped (dt too large, particle hash %u)\n", p->hash);
        return;
    }

    Omega->x *= scale;
    Omega->y *= scale;
    Omega->z *= scale;
}

/* ------------------------------------------------------------------------- */
/* Operator kernel                                                           */
/* ------------------------------------------------------------------------- */
void rebx_magnetic_braking(struct reb_simulation* const sim,
                           struct rebx_operator*   const op,
                           const double                   dt)
{
    struct rebx_extras* const rx = sim->extras;
    const int N = sim->N;

    /* Operator‑level parameters */
    double Msun = 1.0;   /* code‑mass units per solar mass  */
    double Rsun = 1.0;   /* code‑length units per solar radius */
    double year = 1.0;   /* code‑time units per Julian year */

    const double* pM = rebx_get_param(rx, op->ap, "mb_Msun");
    const double* pR = rebx_get_param(rx, op->ap, "mb_Rsun");
    const double* pY = rebx_get_param(rx, op->ap, "mb_year");
    const double* K_ptr = rebx_get_param(rx, op->ap, "mb_K");

    if (pM && isfinite(*pM) && *pM > 0.0) Msun = *pM;
    if (pR && isfinite(*pR) && *pR > 0.0) Rsun = *pR;
    if (pY && isfinite(*pY) && *pY > 0.0) year = *pY;

    const double K_cgs = (K_ptr && isfinite(*K_ptr) && *K_ptr > 0.0)
                            ? *K_ptr : 2.7e47;

    /* Convert prefactor to code units */
    const double Msun_cgs = 1.98847e33;
    const double Rsun_cgs = 6.957e10;
    const double year_cgs = 3.15576e7;
    const double M_unit   = Msun_cgs / Msun;      /* [g/code-mass] */
    const double L_unit   = Rsun_cgs / Rsun;      /* [cm/code-length] */
    const double T_unit   = year_cgs / year;      /* [s/code-time] */
    const double K = K_cgs / (pow(M_unit,1.5) * pow(L_unit,1.5) * T_unit);

    if (dt <= 0.0 || !isfinite(dt)) return;

    double Rossby_sat = 0.1; /* default critical Rossby number */
    const double* Ro_ptr = rebx_get_param(rx, op->ap, "mb_Rossby_sat");
    if (Ro_ptr && isfinite(*Ro_ptr) && *Ro_ptr > 0.0) Rossby_sat = *Ro_ptr;

    for (int i = 0; i < N; i++){
        apply_magnetic_brake(&sim->particles[i], rx, K, dt, Rossby_sat);
    }
}
