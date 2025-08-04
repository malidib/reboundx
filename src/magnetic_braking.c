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
 * I              (double, required)     – moment of inertia
 * Omega          (reb_vec3d, required)  – spin angular‑velocity vector
 *
 * Operator‑level parameters
 * -------------------------
 * mb_K           (double, optional)     – braking constant (default 2.7e47 cgs)
 *
 * NOTE:  K assumes R and M are supplied in the same units used to calibrate K.
 *        If your simulation units differ, rescale mb_K accordingly.
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"

static inline void apply_magnetic_brake(struct reb_particle*    const p,
                                        struct rebx_extras*     const rx,
                                        const double            K,
                                        const double            dt)
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
    const double  omega_sat = (sat_ptr && *sat_ptr > 0.0 && isfinite(*sat_ptr))
                              ? *sat_ptr : INFINITY;

    /* --------------------------- torque magnitude ----------------------- */
    double omega_term = omega*omega*omega;          /* unsaturated default  */
    if (omega > omega_sat)                          /* saturated regime     */
        omega_term = omega * omega_sat * omega_sat;

    const double torque = -K * sqrt(R) * (1.0 / sqrt(M)) * omega_term;  /* cgs‑like */

    /* --------------------------- apply dΩ/dt ---------------------------- */
    const double domega_dt = torque / (*I_ptr);     /* scalar rate (negative) */
    const double scale     = 1.0 + (domega_dt / omega)*dt;

    /* Prevent overshoot that could flip the spin */
    if (scale <= 0.0) return;

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

    const double* K_ptr = rebx_get_param(rx, op->ap, "mb_K");
    const double  K     = (K_ptr && isfinite(*K_ptr) && *K_ptr > 0.0)
                          ? *K_ptr : 2.7e47;      /* default from literature */

    if (dt <= 0.0 || !isfinite(dt)) return;

    for (int i = 0; i < N; i++){
        apply_magnetic_brake(&sim->particles[i], rx, K, dt);
    }
}
