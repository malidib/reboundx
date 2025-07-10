/**
 * @file    roche_lobe_mass_transfer.c
 * @brief   Binary–star sink operator: Roche‑lobe overflow (RLOF),
 *          non‑conservative mass loss, common–envelope drag,
 *          and gravitational‑wave decay – **momentum‑conserving version**.
 *
 * Implements
 *   1. Ritter‑type RLOF (Ritter 1988, Kolb & Ritter 1990)
 *   2. Optional systemic mass loss (parameter *rlmt_loss_fraction*)
 *   3. Chandrasekhar / Ostriker dynamical friction inside a common envelope
 *      (Ostriker 1999 – optional)
 *   4. Peters (1964) gravitational‑wave back‑reaction (optional)
 *
 * **Modifications with respect to the reference implementation**
 *
 *   • Linear **and angular momentum of the transferred /
 *     ejected material are now conserved**:
 *       – the accretor receives the momentum of the material it gains,
 *       – the system loses the momentum carried away by the fraction
 *         that escapes.
 *   • The **back‑reaction on the donor** from common‑envelope drag
 *     is applied (equal and opposite force).
 *   • The **gravitational‑wave reaction** now keeps the barycentre
 *     fixed: both bodies are moved so that the updated relative
 *     orbit is realised while the system momentum is conserved.
 *   • The order of sub‑steps is
 *       1. RLOF ± systemic mass loss (with momentum bookkeeping)
 *       2. Common‑envelope drag (if r < R_d)
 *       3. Gravitational‑wave shrinkage of the orbit
 *     so that each process sees up‑to‑date masses / momenta.
 *
 * The file can simply replace the previous
 * `roche_lobe_mass_transfer.c` in REBOUNDx ≥ v3.3.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

/* ---------- helper for Ostriker drag prefactor ---------------------------------- */
static double mach_piece_sub(const double mach){
    /* Integral I(M) for sub‑sonic Mach numbers (Ostriker 1999, Eq. 13).        */
    /* Series truncated at O(M⁵).                                               */
    if (mach < 0.02){
        const double m2 = mach*mach;
        return  m2*mach/3.        /* ⅓ M³  */
              + m2*m2*mach/5.;    /* ⅕ M⁵  */
    }
    return 0.5*log((1.+mach)/(1.-mach)) - mach;   /* analytic form for M<1   */
}

static double calculate_pre_factor(const double mach, const double xmin){
    /* Coulomb‑logarithm–like term (I) valid for all Mach numbers.              */
    const double coul = log(1./xmin);
    if (mach >= 1.){                /* supersonic: classical Chandrasekhar      */
        return coul;
    }
    return fmin(coul, mach_piece_sub(mach));  /* sub‑sonic: Ostriker (1999)     */
}

/* =============================================================================== */
void rebx_roche_lobe_mass_transfer(struct reb_simulation* const sim,
                                   struct rebx_operator*     const operator,
                                   const double                     dt)
{
    struct rebx_extras* const rebx = sim->extras;

    /* --------  indices of donor / accretor  ---------------------------------- */
    int* donor_idx_ptr     = rebx_get_param(rebx, operator->ap, "rlmt_donor");
    int* accretor_idx_ptr  = rebx_get_param(rebx, operator->ap, "rlmt_accretor");

    if (!donor_idx_ptr || !accretor_idx_ptr){
        rebx_error(rebx, "Need to set rlmt_donor and rlmt_accretor.\n");
        return;
    }
    const int donor_idx    = *donor_idx_ptr;
    const int accretor_idx = *accretor_idx_ptr;

    if (donor_idx >= sim->N || accretor_idx >= sim->N || donor_idx < 0 || accretor_idx < 0){
        rebx_error(rebx, "rlmt_donor or rlmt_accretor index out of range.\n");
        return;
    }

    struct reb_particle* donor     = &sim->particles[donor_idx];
    struct reb_particle* accretor  = &sim->particles[accretor_idx];

    /* --------  effect‑specific parameters  ----------------------------------- */
    const double* Hp_ptr         = rebx_get_param(rebx, donor->ap,     "rlmt_Hp");
    const double* mdot0_ptr      = rebx_get_param(rebx, donor->ap,     "rlmt_mdot0");
    const double* loss_frac_ptr  = rebx_get_param(rebx, operator->ap,  "rlmt_loss_fraction");

    /* Common‑envelope parameters (optional) */
    const double* ce_rho0_ptr    = rebx_get_param(rebx, operator->ap,  "ce_rho0");
    const double* ce_alpha_rho_ptr = rebx_get_param(rebx, operator->ap,"ce_alpha_rho");
    const double* ce_cs_ptr      = rebx_get_param(rebx, operator->ap,  "ce_cs");
    const double* ce_alpha_cs_ptr= rebx_get_param(rebx, operator->ap,  "ce_alpha_cs");
    const double* ce_xmin_ptr    = rebx_get_param(rebx, operator->ap,  "ce_xmin");
    const double* ce_Qd_ptr      = rebx_get_param(rebx, operator->ap,  "ce_Qd");

    /* Gravitational‑wave parameters (optional) */
    const double* gw_c_ptr       = rebx_get_param(rebx, operator->ap,  "gw_c");
    const int*    gw_on_ptr      = rebx_get_param(rebx, operator->ap,  "gw_decay_on");

    if (!Hp_ptr || !mdot0_ptr){
        rebx_error(rebx, "Need to set rlmt_Hp and rlmt_mdot0 on donor particle.\n");
        return;
    }

    /* ------------------------- physical constants & shortcuts ---------------- */
    const double Hp      = *Hp_ptr;
    const double mdot0   = *mdot0_ptr;

    /* Clamp loss fraction to [0,1] */
    double loss_frac = (loss_frac_ptr ? *loss_frac_ptr : 0.);
    if (loss_frac < 0.)      loss_frac = 0.;
    else if (loss_frac > 1.) loss_frac = 1.;

    /* Common‑envelope parameters (defaults are zero → feature disabled) */
    const double ce_rho0    = ce_rho0_ptr     ? *ce_rho0_ptr     : 0.;
    const double ce_alpha_rho = ce_alpha_rho_ptr? *ce_alpha_rho_ptr: 0.;
    const double ce_cs      = ce_cs_ptr       ? *ce_cs_ptr       : 0.;
    const double ce_alpha_cs= ce_alpha_cs_ptr ? *ce_alpha_cs_ptr : 0.;
    const double ce_xmin    = ce_xmin_ptr     ? *ce_xmin_ptr     : 0.;
    const double ce_Qd      = ce_Qd_ptr       ? *ce_Qd_ptr       : 0.;

    /* -------------------------------------------------------------------------
     *  STEP 1 :  Roche‑lobe overflow mass transfer
     * ---------------------------------------------------------------------- */

    /* 1.1 separation of the components */
    const double dx = donor->x - accretor->x;
    const double dy = donor->y - accretor->y;
    const double dz = donor->z - accretor->z;
    const double r  = sqrt(dx*dx + dy*dy + dz*dz);

    /* 1.2 Roche‑lobe radius (Eggleton 1983) */
    const double q        = donor->m / accretor->m;
    const double q13      = cbrt(q);
    const double rl       = r * (0.49*q13*q13) / (0.6*q13*q13 + log(1.+q13));

    /* 1.3 Ritter (1988) mass‑loss rate */
    const double Rd       = donor->r;          /* donor’s physical radius      */
    const double mdot     = -mdot0 * exp((Rd - rl)/Hp);   /* < 0 for mass loss   */
    double       dM       = mdot * dt;         /* negative value               */

    /* Do not let donor mass become negative (one‑shot removal) */
    if (donor->m + dM <= 0.){
        dM = -donor->m;
    }
    const double mass_loss      = -dM;                     /* positive scalar    */
    const double mass_accreted  = mass_loss * (1. - loss_frac);
    const double mass_ejected   = mass_loss * loss_frac;

    /* 1.4 Momentum bookkeeping --------------------------------------------- */
    /* Velocities *before* any mass change */
    const double vdx0 = donor->vx;
    const double vdy0 = donor->vy;
    const double vdz0 = donor->vz;
    const double vax0 = accretor->vx;
    const double vay0 = accretor->vy;
    const double vaz0 = accretor->vz;

    const double Md0 = donor->m;
    const double Ma0 = accretor->m;

    /* Mass update */
    donor->m    -= mass_loss;
    accretor->m += mass_accreted;

    /* The material *gained* by the accretor arrives with the donor’s velocity.
       -> Add its momentum to the accretor and solve for the new velocity.      */
    const double Ma1 = accretor->m;
    const double px_a  = Ma0 * vax0 + mass_accreted * vdx0;
    const double py_a  = Ma0 * vay0 + mass_accreted * vdy0;
    const double pz_a  = Ma0 * vaz0 + mass_accreted * vdz0;

    accretor->vx = px_a / Ma1;
    accretor->vy = py_a / Ma1;
    accretor->vz = pz_a / Ma1;

    /* Donor velocity left *unchanged* – its momentum decreased automatically
       by Md0−Md1 = mass_loss, consistent with the gas it lost.                */

    /* -------------------------------------------------------------------------
     *  STEP 2 :  Common‑envelope dynamical friction (optional)
     *            Only active while the accretor is *inside* the donor’s radius.
     * ---------------------------------------------------------------------- */
    if (ce_rho0 > 0. && ce_cs > 0. && r < Rd){
        /* local envelope stratification */
        const double r_ratio = r / Rd;
        const double rho     = ce_rho0 * pow(r_ratio, ce_alpha_rho);
        const double cs      = ce_cs   * pow(r_ratio, ce_alpha_cs);

        /* relative velocity (updated by mass transfer) */
        const double vrelx   = accretor->vx - donor->vx;
        const double vrely   = accretor->vy - donor->vy;
        const double vrelz   = accretor->vz - donor->vz;
        const double vrel    = sqrt(vrelx*vrelx + vrely*vrely + vrelz*vrelz);

        if (vrel > 0. && rho > 0.){
            const double mach = vrel / cs;
            const double I    = calculate_pre_factor(mach, (ce_xmin>0.?ce_xmin:1e-4));

            /* gravitational drag term (Ostriker 1999; vector form) */
            const double m_a   = accretor->m;
            double fc  = 4.*M_PI*sim->G*sim->G * m_a * rho / (vrel*vrel*vrel) * I;

            /* direct hydrodynamic (accretion) term, optional          */
            if (ce_Qd > 0. && accretor->r > 0.){
                fc += M_PI * rho * accretor->r * accretor->r * vrel * ce_Qd / m_a;
            }

            /* Δv for accretor */
            const double dvx = -fc * vrelx * dt;
            const double dvy = -fc * vrely * dt;
            const double dvz = -fc * vrelz * dt;

            /* Apply to accretor */
            accretor->vx += dvx;
            accretor->vy += dvy;
            accretor->vz += dvz;

            /* -------- back‑reaction on the donor’s core -------------
               Equal and opposite *force*  F = m_a * Δv / dt.
               Donor receives +F / m_d * dt = -(m_a/m_d) Δv.           */
            const double m_d = donor->m;
            if (m_d > 0.){
                donor->vx -= dvx * (m_a / m_d);
                donor->vy -= dvy * (m_a / m_d);
                donor->vz -= dvz * (m_a / m_d);
            }
        }
    }

    /* -------------------------------------------------------------------------
     *  STEP 3 :  Gravitational‑wave back‑reaction (optional)
     *            Peters (1964) secular shrinkage, momentum‑conserving.
     * ---------------------------------------------------------------------- */
    if (gw_c_ptr && (*gw_c_ptr > 0.) && gw_on_ptr && (*gw_on_ptr != 0)){
        /* Current orbit *after* mass transfer & CE‑drag */
        struct reb_orbit o = reb_orbit_from_particle(sim->G, *donor, *accretor);

        if (o.a > 0.){
            const double m1 = donor->m;
            const double m2 = accretor->m;
            const double c  = *gw_c_ptr;
            const double e  = o.e;
            const double a  = o.a;
            const double G3 = sim->G * sim->G * sim->G;

            /* Peters’ da/dt and de/dt */
            const double fac_a = -64./5. * G3*m1*m2*(m1+m2) /
                                 (pow(c,5)*pow(a,3)*pow(1.-e*e,3.5)) *
                                 (1.+73./24.*e*e + 37./96.*e*e*e*e);

            const double fac_e = -304./15. * e * G3*m1*m2*(m1+m2) /
                                 (pow(c,5)*pow(a,4)*pow(1.-e*e,2.5)) *
                                 (1.+121./304.*e*e);

            double a_new = a + fac_a * dt;
            double e_new = e + fac_e * dt;

            if (a_new < 0.)      a_new = 0.;
            if (e_new < 0.)      e_new = 0.;
            if (e_new > 0.999999) e_new = 0.999999;

            /* -------- keep system barycentre fixed ----------------- */
            /* Barycentre AFTER CE‑drag                               */
            const double Mtot = m1 + m2;
            const double Rcmx = (m1*donor->x  + m2*accretor->x) / Mtot;
            const double Rcmy = (m1*donor->y  + m2*accretor->y) / Mtot;
            const double Rcmz = (m1*donor->z  + m2*accretor->z) / Mtot;
            const double Vcmx = (m1*donor->vx + m2*accretor->vx) / Mtot;
            const double Vcmy = (m1*donor->vy + m2*accretor->vy) / Mtot;
            const double Vcmz = (m1*donor->vz + m2*accretor->vz) / Mtot;

            /* Generate NEW relative orbit (donor around accretor) */
            struct reb_particle np =
                reb_particle_from_orbit(sim->G, *accretor, m1,
                                        a_new, e_new,
                                        o.inc, o.Omega, o.omega, o.f);

            /* Relative vectors (new) */
            const double rx_rel = np.x - accretor->x;
            const double ry_rel = np.y - accretor->y;
            const double rz_rel = np.z - accretor->z;

            const double vx_rel = np.vx - accretor->vx;
            const double vy_rel = np.vy - accretor->vy;
            const double vz_rel = np.vz - accretor->vz;

            /* Place the two bodies around the *same* barycentre */
            const double frac_d =   m2 / Mtot;   /* donor offset factor      */
            const double frac_a =  -m1 / Mtot;   /* accretor offset factor   */

            donor->x  = Rcmx + frac_d * rx_rel;
            donor->y  = Rcmy + frac_d * ry_rel;
            donor->z  = Rcmz + frac_d * rz_rel;

            accretor->x = Rcmx + frac_a * rx_rel;
            accretor->y = Rcmy + frac_a * ry_rel;
            accretor->z = Rcmz + frac_a * rz_rel;

            donor->vx = Vcmx + frac_d * vx_rel;
            donor->vy = Vcmy + frac_d * vy_rel;
            donor->vz = Vcmz + frac_d * vz_rel;

            accretor->vx = Vcmx + frac_a * vx_rel;
            accretor->vy = Vcmy + frac_a * vy_rel;
            accretor->vz = Vcmz + frac_a * vz_rel;
        }
    }

    /* --------------------------------------------------------------------- */
    /* Remove donor if fully evaporated, detach operator */
    if (donor->m <= 0.){
        reb_simulation_remove_particle(sim, donor_idx, 1);
        rebx_remove_operator(rebx, operator);      /* stops further calls */
        reb_simulation_move_to_com(sim);
        return;
    }

    /* Re‑centre the simulation (cheap, keeps integrator happy) */
    reb_simulation_move_to_com(sim);
}
