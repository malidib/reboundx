/* ============================================================================
 * @file    roche_lobe_mass_transfer.c
 * @brief   Roche‑Lobe Overflow (RLOF) + Common‑Envelope (CE) operator
 *          with momentum‑exact mass transfer, controlled j‑loss,
 *          optional CE reaction, adaptive sub‑stepping, and diagnostics.
 *
 *          2025‑08‑15  (fully corrected)
 *
 * Operator parameters (all optional unless stated)
 * ------------------------------------------------
 * rlmt_donor            (double, required)  – donor particle index (0‑based)
 * rlmt_accretor         (double, required)  – accretor particle index (0‑based)
 * rlmt_loss_fraction    (double)            – wind fraction f_loss in [0,1]     (default 0)
 * jloss_mode            (double)            – 0 donor wind (v_loss = v_d)
 *                                            1 isotropic re‑emission (v_a)
 *                                            2 Jeans (v_CM)
 *                                            3 scale × j_orb (target j)         (default 0)
 * jloss_factor          (double)            – scale for mode 3 (default 1.0)
 * rlmt_skip_in_CE       (double, bool)      – skip RLOF if r < R_d               (default 1)
 * rlmt_substep_max_dm   (double)            – max |ΔM|/M per sub‑step            (default 1e‑3)
 * rlmt_substep_max_dr   (double)            – max |Δr|/r per sub‑step            (default 5e‑3)
 * rlmt_min_substeps     (double)            – enforce ≥N sub‑steps               (default 3)
 * ce_profile_file       (string, optional)  – ASCII table (s, rho, cs), overrides power‑law
 * ce_kick_cfl           (double)            – |Δv| ≤ cfl × c_s per sub‑step      (default 1.0)
 * ce_reaction_on_donor  (double, bool)      – apply opposite CE kick to donor    (default 0)
 * merge_eps             (double)            – explicit merge radius; default 0.5×min(Rs) with floor
 *
 * Particle parameters (donor unless noted)
 * ---------------------------------------
 * rlmt_Hp               (double, donor)     – pressure scale height H_P
 * rlmt_mdot0            (double, donor)     – reference mass‑loss rate \dot M_0  (>0)
 * rlmt_R_slope          (double, donor)     – d ln R / d ln M exponent α_R       (default 0)
 * rlmt_R_ref_mass       (double, donor)     – reference mass for R(M)            (default donor’s current)
 * rlmt_R_ref_radius     (double, donor)     – reference radius for R(M)          (default donor’s current)
 *
 * CE power‑law (operator scope; used if no table)
 * -----------------------------------------------
 * ce_rho0               (double)            – density normalization
 * ce_alpha_rho          (double)            – density slope
 * ce_cs                 (double)            – sound‑speed normalization
 * ce_alpha_cs           (double)            – sound‑speed slope
 * ce_xmin               (double)            – Coulomb cutoff x_min               (default 1e‑4)
 * ce_Qd                 (double)            – geometric drag coefficient          (default 0)
 *
 * Output diagnostic
 * -----------------
 * rlmt_last_dE          (double, op‑attr)   – total ΔE from last call
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

#ifndef RLMT_EXP_CLAMP
#define RLMT_EXP_CLAMP  80.0     /* prevents exp() overflow in Ritter law        */
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define MAX2(a,b) (( (a) > (b) ) ? (a) : (b))
#define MIN2(a,b) (( (a) < (b) ) ? (a) : (b))

/* --- Optional string param getter (define REBX_ENABLE_STRING_PARAMS if available) --- */
#ifdef REBX_ENABLE_STRING_PARAMS
/* REBOUNDx provides this in newer versions; forward declare to avoid header mismatch. */
extern const char* rebx_get_param_str(struct rebx_extras* rx, struct rebx_param* ap, const char* name);
#endif

/* -------------------------- CE profile table (global) -------------------------- */
/* Kept global for simplicity; loaded at first use if ce_profile_file is set.     */
struct ce_profile {
    int     n;
    double* s;      /* r/R_d */
    double* rho;    /* density */
    double* cs;     /* sound speed */
};
static struct ce_profile ce_tab = {0, NULL, NULL};

/* Load (s, rho, cs) ASCII table. Returns 1 on success, 0 otherwise. */
static int ce_profile_load(const char* fname){
    if(!fname) return 0;
    FILE* f = fopen(fname, "r");
    if(!f){
        fprintf(stderr, "[rlmt] Cannot open CE profile file: %s\n", fname);
        return 0;
    }
    /* free old */
    free(ce_tab.s); free(ce_tab.rho); free(ce_tab.cs);
    ce_tab.n = 0; ce_tab.s = ce_tab.rho = ce_tab.cs = NULL;

    int cap = 1024;
    ce_tab.s   = (double*)malloc((size_t)cap*sizeof(double));
    ce_tab.rho = (double*)malloc((size_t)cap*sizeof(double));
    ce_tab.cs  = (double*)malloc((size_t)cap*sizeof(double));
    if(!ce_tab.s || !ce_tab.rho || !ce_tab.cs){
        fprintf(stderr, "[rlmt] Out of memory reading CE profile.\n");
        fclose(f);
        return 0;
    }

    while(1){
        double x, r, c;
        int nread = fscanf(f, "%lf %lf %lf", &x, &r, &c);
        if(nread != 3) break;
        if(x <= 0.0 || r <= 0.0 || c <= 0.0) continue;  /* require positive for log‑interp */
        if(ce_tab.n == cap){
            cap *= 2;
            double* ns  = (double*)realloc(ce_tab.s,   (size_t)cap*sizeof(double));
            double* nr  = (double*)realloc(ce_tab.rho, (size_t)cap*sizeof(double));
            double* ncs = (double*)realloc(ce_tab.cs,  (size_t)cap*sizeof(double));
            if(!ns || !nr || !ncs){
                fprintf(stderr, "[rlmt] Out of memory expanding CE profile.\n");
                fclose(f);
                return 0;
            }
            ce_tab.s = ns; ce_tab.rho = nr; ce_tab.cs = ncs;
        }
        ce_tab.s  [ce_tab.n] = x;
        ce_tab.rho[ce_tab.n] = r;
        ce_tab.cs [ce_tab.n] = c;
        ce_tab.n++;
    }
    fclose(f);

    if(ce_tab.n < 2){
        fprintf(stderr, "[rlmt] CE profile file has fewer than 2 valid rows; ignoring.\n");
        free(ce_tab.s); free(ce_tab.rho); free(ce_tab.cs);
        ce_tab.n = 0; ce_tab.s = ce_tab.rho = ce_tab.cs = NULL;
        return 0;
    }
    return 1;
}

/* Log‑log linear interpolation (positive‑definite). */
static inline void ce_interp(double s, double* rho_out, double* cs_out){
    if(!rho_out || !cs_out) return;
    if(ce_tab.n < 2){ *rho_out = 0.0; *cs_out = 0.0; return; }
    if(s <= ce_tab.s[0]){ *rho_out = ce_tab.rho[0]; *cs_out = ce_tab.cs[0]; return; }
    if(s >= ce_tab.s[ce_tab.n-1]){ *rho_out = ce_tab.rho[ce_tab.n-1]; *cs_out = ce_tab.cs[ce_tab.n-1]; return; }

    /* binary search */
    int i = 0, j = ce_tab.n - 1;
    while(j - i > 1){
        int m = (i + j) >> 1;
        if(s < ce_tab.s[m]) j = m; else i = m;
    }
    double t = (s - ce_tab.s[i]) / (ce_tab.s[j] - ce_tab.s[i]);
    /* positive inputs guaranteed */
    double lrho = log(ce_tab.rho[i])*(1.0 - t) + log(ce_tab.rho[j])*t;
    double lcs  = log(ce_tab.cs [i])*(1.0 - t) + log(ce_tab.cs [j])*t;
    *rho_out = exp(lrho);
    *cs_out  = exp(lcs);
}

/* ------------------------------ small helpers ------------------------------ */
static inline void cross3(const double ax, const double ay, const double az,
                          const double bx, const double by, const double bz,
                          double* rx, double* ry, double* rz){
    *rx = ay*bz - az*by;
    *ry = az*bx - ax*bz;
    *rz = ax*by - ay*bx;
}
static inline double dot3(const double ax, const double ay, const double az,
                          const double bx, const double by, const double bz){
    return ax*bx + ay*by + az*bz;
}

/* Build any unit vector perpendicular to n = (nx,ny,nz). */
static void unit_perp_to(const double nx, const double ny, const double nz,
                         double* ex, double* ey, double* ez){
    /* try cross with z‑axis first, then x‑axis */
    double cx = ny*1.0 - nz*0.0;  /* n × ez */
    double cy = nz*0.0 - nx*1.0;
    double cz = nx*0.0 - ny*0.0;
    double nrm = sqrt(cx*cx + cy*cy + cz*cz);
    if(nrm < 1e-12){
        cx = ny*0.0 - nz*0.0;     /* n × ex */
        cy = nz*1.0 - nx*0.0;
        cz = nx*0.0 - ny*1.0;
        nrm = sqrt(cx*cx + cy*cy + cz*cz);
        if(nrm < 1e-12){ *ex = 1.0; *ey = 0.0; *ez = 0.0; return; }
    }
    *ex = cx/nrm; *ey = cy/nrm; *ez = cz/nrm;
}

/* Ostriker low‑Mach series piece; analytic for M<1 else handled in I_prefactor. */
static double mach_piece_sub(const double M){
    if(M < 0.02){
        const double m2 = M*M;
        return m2*M/3. + m2*m2*M/5.;       /* ⅓ M^3 + ⅕ M^5 */
    }
    return 0.5*log((1.+M)/(1.-M)) - M;     /* analytic for 0.02 ≤ M < 1 */
}
static double I_prefactor(const double M, const double xmin){
    const double coul = log(1.0/xmin);
    if(M >= 1.0) return coul;
    /* cap by Coulomb log to keep continuous across regimes */
    return MIN2(coul, mach_piece_sub(M));
}

/* ========================================================================= */
void rebx_roche_lobe_mass_transfer(struct reb_simulation* const sim,
                                   struct rebx_operator*     const op,
                                   const double                     dt_req)
{
    if(!(isfinite(dt_req) && dt_req > 0.0)) return;

    struct rebx_extras* const rx = sim->extras;

    /* ---------------------- required indices (as doubles) ---------------------- */
    const double* d_donor = rebx_get_param(rx, op->ap, "rlmt_donor");
    const double* d_accr  = rebx_get_param(rx, op->ap, "rlmt_accretor");
    if(!d_donor || !d_accr){
        reb_simulation_error(sim, "[rlmt] Need rlmt_donor and rlmt_accretor (double indices).");
        return;
    }
    const int donor_idx = (int) llround(*d_donor);
    const int acc_idx   = (int) llround(*d_accr);
    if(donor_idx < 0 || acc_idx < 0 || donor_idx >= sim->N || acc_idx >= sim->N || donor_idx == acc_idx){
        reb_simulation_error(sim, "[rlmt] rlmt_donor/rlmt_accretor out of range or identical.");
        return;
    }
    struct reb_particle* d = &sim->particles[donor_idx];
    struct reb_particle* a = &sim->particles[acc_idx];

    if(!(d->m > 0.0 && a->m > 0.0)) return;  /* nothing to do if any zero/neg mass */

    /* --------------------------- operator parameters --------------------------- */
    const double* p_loss   = rebx_get_param(rx, op->ap, "rlmt_loss_fraction");
    double f_loss = p_loss ? *p_loss : 0.0;
    if(!isfinite(f_loss)) f_loss = 0.0;
    if(f_loss < 0.0) f_loss = 0.0; if(f_loss > 1.0) f_loss = 1.0;

    const double* p_jmode  = rebx_get_param(rx, op->ap, "jloss_mode");
    const double* p_jfac   = rebx_get_param(rx, op->ap, "jloss_factor");
    const int     jloss_mode   = p_jmode ? (int) llround(*p_jmode) : 0;
    const double  jloss_factor = p_jfac  ? *p_jfac : 1.0;

    const double* p_skipCE = rebx_get_param(rx, op->ap, "rlmt_skip_in_CE");
    const int     skip_in_CE = p_skipCE ? (int) llround(*p_skipCE) : 1;

    const double* p_dmmax  = rebx_get_param(rx, op->ap, "rlmt_substep_max_dm");
    const double* p_drmax  = rebx_get_param(rx, op->ap, "rlmt_substep_max_dr");
    const double* p_nmin   = rebx_get_param(rx, op->ap, "rlmt_min_substeps");
    const double* p_cfl    = rebx_get_param(rx, op->ap, "ce_kick_cfl");
    const double* p_merge  = rebx_get_param(rx, op->ap, "merge_eps");
    const double* p_react  = rebx_get_param(rx, op->ap, "ce_reaction_on_donor");

    const double dm_max_frac = p_dmmax ? *p_dmmax : 1e-3;
    const double dr_max_frac = p_drmax ? *p_drmax : 5e-3;
    const int    min_steps   = p_nmin  ? (int) llround(*p_nmin) : 3;
    const double kick_cfl    = p_cfl   ? *p_cfl : 1.0;
    const int    ce_react_on_donor = p_react ? (int) llround(*p_react) : 0;

    /* optional CE profile table (loaded once) */
    static int ce_loaded = 0;
#ifdef REBX_ENABLE_STRING_PARAMS
    if(!ce_loaded){
        const char* ce_file = rebx_get_param_str(rx, op->ap, "ce_profile_file");
        if(ce_file && ce_file[0] != '\0'){
            ce_loaded = ce_profile_load(ce_file);
            if(!ce_loaded){
                fprintf(stderr, "[rlmt] CE table load failed; falling back to power‑law CE.\n");
            }
        } else {
            ce_loaded = -1; /* tried; no file */
        }
    }
#endif

    /* --------------------------- initial merge guard --------------------------- */
    const double dx0 = d->x - a->x;
    const double dy0 = d->y - a->y;
    const double dz0 = d->z - a->z;
    const double r2_0 = dx0*dx0 + dy0*dy0 + dz0*dz0;
    const double r0   = sqrt(r2_0);

    double merge_eps = p_merge ? *p_merge : 0.0;
    if(!(merge_eps > 0.0)){
        double eps_default = 0.0;
        if(d->r > 0.0 && a->r > 0.0) eps_default = 0.5 * MIN2(d->r, a->r);
        /* positive floor, e.g., 1e-6 of current separation if radii lack info */
        if(!(eps_default > 0.0)) eps_default = 1e-6 * MAX2(r0, 1.0);
        merge_eps = eps_default;
    }

    if(!(isfinite(r2_0)) || r2_0 <= merge_eps*merge_eps){
        /* Merge immediately */
        const double msum = d->m + a->m;
        if(msum > 0.0){
            d->x  = (d->m*d->x  + a->m*a->x ) / msum;
            d->y  = (d->m*d->y  + a->m*a->y ) / msum;
            d->z  = (d->m*d->z  + a->m*a->z ) / msum;
            d->vx = (d->m*d->vx + a->m*a->vx) / msum;
            d->vy = (d->m*d->vy + a->m*a->vy) / msum;
            d->vz = (d->m*d->vz + a->m*a->vz) / msum;
            d->m  = msum;
            d->r  = MAX2(d->r, a->r);
        } else {
            d->m = 0.0;
        }
        reb_simulation_remove_particle(sim, acc_idx, 1);
        if(sim->N < 2) rebx_remove_operator(rx, op);
        reb_simulation_move_to_com(sim);
        return;
    }

    /* --------------------------- sub‑stepping loop ---------------------------- */
    double t_left = dt_req;
    int    steps  = 0;
    double dE_total = 0.0;  /* diagnostic only */

    while(t_left > 0.0 && (steps < min_steps || t_left > 1e-14*dt_req)){
        double dt = t_left / (double)((min_steps - steps) > 0 ? (min_steps - steps) : 1);

        /* Separation & relative kinematics */
        double dx = d->x - a->x;
        double dy = d->y - a->y;
        double dz = d->z - a->z;
        double r  = sqrt(dx*dx + dy*dy + dz*dz);
        if(!(r > 0.0)) break;

        const double nx = dx / r, ny = dy / r, nz = dz / r;

        const double vrelx = a->vx - d->vx;
        const double vrely = a->vy - d->vy;
        const double vrelz = a->vz - d->vz;
        const double vrel  = sqrt(vrelx*vrelx + vrely*vrely + vrelz*vrelz);

        /* donor RLOF parameters (required) */
        const double* Hp_ptr    = rebx_get_param(rx, d->ap, "rlmt_Hp");
        const double* mdot0_ptr = rebx_get_param(rx, d->ap, "rlmt_mdot0");
        if(!Hp_ptr || !mdot0_ptr){
            reb_simulation_error(sim, "[rlmt] Donor needs rlmt_Hp and rlmt_mdot0.");
            return;
        }
        const double Hp    = *Hp_ptr;
        const double mdot0 = *mdot0_ptr;
        if(!(Hp > 0.0 && mdot0 > 0.0)){
            reb_simulation_error(sim, "[rlmt] rlmt_Hp and rlmt_mdot0 must be positive.");
            return;
        }

        /* Eggleton Roche radius (donor as primary) */
        const double q    = d->m / a->m;
        const double q13  = cbrt(q);
        const double q23  = q13*q13;
        const double RL   = r * (0.49*q23) / (0.6*q23 + log(1.0 + q13));

        /* Ritter mass‑loss estimate for step limiting */
        double expo = (d->r - RL) / Hp;
        if(expo > RLMT_EXP_CLAMP) expo = RLMT_EXP_CLAMP;
        double mdot_est = -mdot0 * exp(expo);  /* <0 when overflowing */

        if(mdot_est != 0.0){
            double dt_lim = dm_max_frac * MAX2(d->m, 1e-30) / fabs(mdot_est);
            if(dt > dt_lim) dt = dt_lim;
        }
        /* safety dt based on relative motion */
        double dt_vel = dr_max_frac * r / MAX2(vrel, 1e-15);
        if(dt > dt_vel) dt = dt_vel;

        if(dt > t_left) dt = t_left;

        /* Energy before step (diagnostic only) */
        const double Md0 = d->m, Ma0 = a->m;
        const double E_before =
            0.5*Md0*(d->vx*d->vx + d->vy*d->vy + d->vz*d->vz) +
            0.5*Ma0*(a->vx*a->vx + a->vy*a->vy + a->vz*a->vz) -
            sim->G*Md0*Ma0/r;

        /* ===================================================================== */
        /* STEP 1 – Roche‑lobe overflow (skip if inside CE and skip_in_CE==1)    */
        /* ===================================================================== */
        if(!(skip_in_CE && r < d->r)){
            /* Conservative mass loss from donor */
            expo = (d->r - RL) / Hp; if(expo > RLMT_EXP_CLAMP) expo = RLMT_EXP_CLAMP;
            const double mdot = -mdot0 * exp(expo);   /* <0 => donor losing mass */
            double dM   = mdot * dt;                  /* negative */
            if(d->m + dM <= 0.0) dM = -d->m + (-1e-30);  /* leave tiny positive */

            const double m_loss = -dM;                /* >0: donor mass decrease */
            const double m_wind = f_loss * m_loss;    /* escapes system */
            const double m_acc  = m_loss - m_wind;    /* accreted internally */

            /* Compute v_loss according to jloss_mode */
            double vx_loss, vy_loss, vz_loss;
            if(jloss_mode == 1){
                vx_loss = a->vx; vy_loss = a->vy; vz_loss = a->vz;
            } else if(jloss_mode == 2){
                const double Mtot = Md0 + Ma0;
                vx_loss = (Md0*d->vx + Ma0*a->vx) / Mtot;
                vy_loss = (Md0*d->vy + Ma0*a->vy) / Mtot;
                vz_loss = (Md0*d->vz + Ma0*a->vz) / Mtot;
            } else if(jloss_mode == 3){
                /* choose tangential unit vector and set |Δv| to achieve target j */
                double exu, eyu, ezu; unit_perp_to(nx, ny, nz, &exu, &eyu, &ezu);
                /* j_orb = |r × v| / μ  with r = donor - accretor, v = v_rel */
                double Lx, Ly, Lz; cross3(dx, dy, dz, vrelx, vrely, vrelz, &Lx, &Ly, &Lz);
                const double Lmag = sqrt(Lx*Lx + Ly*Ly + Lz*Lz) + 1e-99;
                const double mu   = (Md0*Ma0) / (Md0 + Ma0);
                const double j_orb = Lmag / mu;
                const double j_target = jloss_factor * j_orb;
                const double fac = j_target / r;  /* speed to achieve j_target at radius r (donor) */
                vx_loss = d->vx + fac*exu;
                vy_loss = d->vy + fac*eyu;
                vz_loss = d->vz + fac*ezu;
            } else {
                /* mode 0: donor wind */
                vx_loss = d->vx; vy_loss = d->vy; vz_loss = d->vz;
            }

            /* --- Mass updates --- */
            const double Md1 = Md0 - m_loss;       /* donor mass after loss */
            const double Ma1 = Ma0 + m_acc;        /* accretor mass after accretion */
            d->m = Md1; a->m = Ma1;

            /* --- Internal conservative accretion: momentum‑exact update --- */
            if(m_acc > 0.0){
                /* Add m_acc at donor velocity to accretor */
                a->vx = (Ma0*a->vx + m_acc*d->vx) / Ma1;
                a->vy = (Ma0*a->vy + m_acc*d->vy) / Ma1;
                a->vz = (Ma0*a->vz + m_acc*d->vz) / Ma1;
                /* Donor velocity unchanged; donor momentum reduced implicitly by mass loss. */
            }

            /* --- External wind: remove momentum −m_wind v_loss via uniform COM shift --- */
            if(m_wind > 0.0){
                const double Mtot1 = Md1 + Ma1;
                const double dVx = -(m_wind * vx_loss) / Mtot1;
                const double dVy = -(m_wind * vy_loss) / Mtot1;
                const double dVz = -(m_wind * vz_loss) / Mtot1;
                d->vx += dVx; d->vy += dVy; d->vz += dVz;
                a->vx += dVx; a->vy += dVy; a->vz += dVz;
            }

            /* --- Conservative ΔL correction for transferred (accreted) mass --- */
            if(m_acc > 0.0){
                /* R_cm after mass update */
                const double Rcmx = (Md1*d->x + Ma1*a->x) / (Md1 + Ma1);
                const double Rcmy = (Md1*d->y + Ma1*a->y) / (Md1 + Ma1);
                const double Rcmz = (Md1*d->z + Ma1*a->z) / (Md1 + Ma1);
                const double rdx = d->x - Rcmx, rdy = d->y - Rcmy, rdz = d->z - Rcmz;
                const double rax = a->x - Rcmx, ray = a->y - Rcmy, raz = a->z - Rcmz;

                /* ΔL_acc = m_acc * [(r_a − r_d) × v_d]; we apply −ΔL_acc as a pure torque */
                double rdXvd_x, rdXvd_y, rdXvd_z; cross3(rdx, rdy, rdz, d->vx, d->vy, d->vz, &rdXvd_x, &rdXvd_y, &rdXvd_z);
                double raXvd_x, raXvd_y, raXvd_z; cross3(rax, ray, raz, d->vx, d->vy, d->vz, &raXvd_x, &raXvd_y, &raXvd_z);
                double dLx = m_acc * (raXvd_x - rdXvd_x);
                double dLy = m_acc * (raXvd_y - rdXvd_y);
                double dLz = m_acc * (raXvd_z - rdXvd_z);

                /* Apply δv_a such that Ma1 (r_a × δv_a) = −ΔL_acc; donor gets opposite to keep ΔP=0 */
                const double ra2 = rax*rax + ray*ray + raz*raz;
                if(ra2 > 0.0 && Ma1 > 0.0 && Md1 > 0.0){
                    /* δv = ( -ΔL × r_a ) / (M_a |r_a|^2) */
                    const double tx = -(dLy*raz - dLz*ray) / (Ma1 * ra2);
                    const double ty = -(dLz*rax - dLx*raz) / (Ma1 * ra2);
                    const double tz = -(dLx*ray - dLy*rax) / (Ma1 * ra2);
                    a->vx += tx; a->vy += ty; a->vz += tz;
                    const double scale = Ma1 / Md1;
                    d->vx -= scale*tx; d->vy -= scale*ty; d->vz -= scale*tz;
                }
            }

            /* --- Wind angular‑momentum removal: enforce ΔL_wind explicitly --- */
            if(m_wind > 0.0){
                /* R_cm after updates */
                const double Rcmx = (d->m*d->x + a->m*a->x) / (d->m + a->m);
                const double Rcmy = (d->m*d->y + a->m*a->y) / (d->m + a->m);
                const double Rcmz = (d->m*d->z + a->m*a->z) / (d->m + a->m);
                double rex, rey, rez;  /* emission point relative to COM */
                if(jloss_mode == 1){      /* accretor wind */
                    rex = a->x - Rcmx; rey = a->y - Rcmy; rez = a->z - Rcmz;
                } else if(jloss_mode == 2){ /* Jeans (COM) */
                    rex = rey = rez = 0.0;
                } else {                   /* donor wind or mode 3 baseline from donor */
                    rex = d->x - Rcmx; rey = d->y - Rcmy; rez = d->z - Rcmz;
                }

                double DLx, DLy, DLz;
                if(jloss_mode == 3){
                    /* Target |ΔL| = m_wind * j_target along orbital L‑hat */
                    double Lx, Ly, Lz; cross3(dx, dy, dz, vrelx, vrely, vrelz, &Lx, &Ly, &Lz);
                    const double Lmag = sqrt(Lx*Lx + Ly*Ly + Lz*Lz) + 1e-99;
                    const double mu   = (d->m*a->m) / (d->m + a->m);
                    const double j_orb = Lmag / mu;
                    const double j_target = jloss_factor * j_orb;
                    const double scaleL = m_wind * j_target / Lmag;
                    DLx = scaleL * Lx; DLy = scaleL * Ly; DLz = scaleL * Lz;
                } else {
                    /* Use ΔL_wind = m_wind * (r_emit × v_loss) */
                    cross3(rex, rey, rez, vx_loss, vy_loss, vz_loss, &DLx, &DLy, &DLz);
                    DLx *= m_wind; DLy *= m_wind; DLz *= m_wind;
                }

                /* Apply −ΔL_wind as a pure torque on the pair */
                const double rax = a->x - Rcmx, ray = a->y - Rcmy, raz = a->z - Rcmz;
                const double ra2 = rax*rax + ray*ray + raz*raz;
                if(ra2 > 0.0 && a->m > 0.0 && d->m > 0.0){
                    const double tx = -(DLy*raz - DLz*ray) / (a->m * ra2);
                    const double ty = -(DLz*rax - DLx*raz) / (a->m * ra2);
                    const double tz = -(DLx*ray - DLy*rax) / (a->m * ra2);
                    a->vx += tx; a->vy += ty; a->vz += tz;
                    const double scale = a->m / d->m;
                    d->vx -= scale*tx; d->vy -= scale*ty; d->vz -= scale*tz;
                }
            }

            /* --- donor mass‑radius relation (optional) --- */
            const double* p_Rslope = rebx_get_param(rx, d->ap, "rlmt_R_slope");
            if(p_Rslope && *p_Rslope != 0.0 && d->m > 0.0){
                const double alpha = *p_Rslope;
                const double* p_Mref = rebx_get_param(rx, d->ap, "rlmt_R_ref_mass");
                const double* p_Rref = rebx_get_param(rx, d->ap, "rlmt_R_ref_radius");
                const double Mref = (p_Mref && *p_Mref > 0.0) ? *p_Mref : Md0;
                const double Rref = (p_Rref && *p_Rref > 0.0) ? *p_Rref : d->r;
                if(Mref > 0.0 && Rref > 0.0){
                    d->r = Rref * pow(d->m / Mref, alpha);
                }
            }
        } /* end RLOF */

        /* ===================================================================== */
        /* STEP 2 – Common‑Envelope dynamical friction (optional)                */
        /* ===================================================================== */
        {
            const double* rho0_ptr = rebx_get_param(rx, op->ap, "ce_rho0");
            const double* cs0_ptr  = rebx_get_param(rx, op->ap, "ce_cs");
            if(r < d->r && ( (rho0_ptr && cs0_ptr) || (ce_tab.n >= 2) )){
                const double* arho_ptr = rebx_get_param(rx, op->ap, "ce_alpha_rho");
                const double* acs_ptr  = rebx_get_param(rx, op->ap, "ce_alpha_cs");
                const double* xmin_ptr = rebx_get_param(rx, op->ap, "ce_xmin");
                const double* Qd_ptr   = rebx_get_param(rx, op->ap, "ce_Qd");

                const double s = r / MAX2(d->r, 1e-30);
                double rho, cs;
                if(ce_tab.n >= 2){
                    ce_interp(s, &rho, &cs);
                } else {
                    const double rho0 = *rho0_ptr;
                    const double cs0  = *cs0_ptr;
                    const double arho = arho_ptr ? *arho_ptr : 0.0;
                    const double acs  = acs_ptr  ? *acs_ptr  : 0.0;
                    rho = rho0 * pow(s, arho);
                    cs  = cs0  * pow(s, acs);
                }
                if(rho > 0.0 && cs > 0.0 && a->m > 0.0){
                    const double xmin = xmin_ptr ? *xmin_ptr : 1e-4;
                    double vrel_eff = MAX2(vrel, 1e-3*cs);
                    const double I = I_prefactor(vrel_eff / cs, xmin);

                    /* Ostriker drag on accretor */
                    const double fc = 4.0 * M_PI * sim->G*sim->G * a->m * rho / (vrel_eff*vrel_eff*vrel_eff) * I;

                    double dvx = -fc * vrelx * dt;
                    double dvy = -fc * vrely * dt;
                    double dvz = -fc * vrelz * dt;

                    /* Optional geometric term */
                    const double* Qd2_ptr = Qd_ptr;
                    if(Qd2_ptr && *Qd2_ptr > 0.0 && a->r > 0.0){
                        const double fc_geom = M_PI * rho * a->r * a->r * vrel / a->m * (*Qd2_ptr);
                        dvx += -fc_geom * vrelx * dt;
                        dvy += -fc_geom * vrely * dt;
                        dvz += -fc_geom * vrelz * dt;
                    }

                    /* CFL‑like limiter */
                    const double dv = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
                    const double cap = kick_cfl * cs;
                    if(dv > cap && dv > 0.0){
                        const double f = cap / dv;
                        dvx *= f; dvy *= f; dvz *= f;
                    }

                    /* Apply to accretor; optional reaction on donor */
                    a->vx += dvx; a->vy += dvy; a->vz += dvz;
                    if(ce_react_on_donor && d->m > 0.0){
                        const double scale = a->m / d->m;
                        d->vx -= scale*dvx; d->vy -= scale*dvy; d->vz -= scale*dvz;
                    }
                }
            }
        }

        /* ===================================================================== */
        /* STEP 3 – Secondary merge guard                                        */
        /* ===================================================================== */
        dx = d->x - a->x; dy = d->y - a->y; dz = d->z - a->z;
        if(dx*dx + dy*dy + dz*dz <= merge_eps*merge_eps){
            const double msum = d->m + a->m;
            if(msum > 0.0){
                d->x  = (d->m*d->x  + a->m*a->x ) / msum;
                d->y  = (d->m*d->y  + a->m*a->y ) / msum;
                d->z  = (d->m*d->z  + a->m*a->z ) / msum;
                d->vx = (d->m*d->vx + a->m*a->vx) / msum;
                d->vy = (d->m*d->vy + a->m*a->vy) / msum;
                d->vz = (d->m*d->vz + a->m*a->vz) / msum;
                d->m  = msum;
                d->r  = MAX2(d->r, a->r);
            } else {
                d->m = 0.0;
            }
            reb_simulation_remove_particle(sim, acc_idx, 1);
            if(sim->N < 2) rebx_remove_operator(rx, op);
            reb_simulation_move_to_com(sim);
            return;
        }

        /* ===================================================================== */
        /* STEP 4 – Purge zero‑mass particles                                    */
        /* ===================================================================== */
        for(int i = sim->N - 1; i >= 0; --i){
            if(sim->particles[i].m <= 0.0){
                reb_simulation_remove_particle(sim, i, 1);
                if(i == donor_idx || i == acc_idx){
                    rebx_remove_operator(rx, op);
                    reb_simulation_move_to_com(sim);
                    return;
                }
            }
        }

        /* Energy after step (diagnostic only) */
        const double E_after =
            0.5*d->m*(d->vx*d->vx + d->vy*d->vy + d->vz*d->vz) +
            0.5*a->m*(a->vx*a->vx + a->vy*a->vy + a->vz*a->vz) -
            sim->G*d->m*a->m / sqrt( (d->x - a->x)*(d->x - a->x) +
                                      (d->y - a->y)*(d->y - a->y) +
                                      (d->z - a->z)*(d->z - a->z) );
        dE_total += (E_after - E_before);

        t_left -= dt;
        steps++;
    }

    /* store diagnostic */
    rebx_set_param_double(rx, &op->ap, "rlmt_last_dE", dE_total);

    /* recentre */
    reb_simulation_move_to_com(sim);
}
