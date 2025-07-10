/**
 * @file    roche_lobe_mass_transfer.c
 * @brief   Binary‑star sink operator: Roche‑lobe overflow, systemic
 *          mass loss, common‑envelope drag, gravitational‑wave decay,
 *          plus automatic object removal & merging.
 *
 * Last update: 2025‑07‑10
 *
 * New in this revision
 *   • Any particle whose mass becomes ≤ 0 is removed immediately.
 *   • Donor and accretor are **merged** when their separation drops
 *     below *merge_eps* (parameter, default = 1 cm).
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

/* ---------- vector helpers ------------------------------------------------- */
static inline void cross3(const double ax, const double ay, const double az,
                          const double bx, const double by, const double bz,
                          double* const rx, double* const ry, double* const rz){
    *rx = ay*bz - az*by;
    *ry = az*bx - ax*bz;
    *rz = ax*by - ay*bx;
}

/* ---------- Ostriker drag helper ------------------------------------------- */
static double mach_piece_sub(const double M){
    if (M < 0.02){
        const double m2 = M*M;
        return  m2*M/3. + m2*m2*M/5.;                /* ⅓ M³ + ⅕ M⁵          */
    }
    return 0.5*log((1.+M)/(1.-M)) - M;               /* analytic M<1         */
}
static double I_prefactor(const double M, const double xmin){
    const double coul = log(1./xmin);
    if (M > 1.0)   return coul;                      /* supersonic           */
    return fmin(coul, mach_piece_sub(M));            /* sub‑sonic            */
}

/* ========================================================================== */
void rebx_roche_lobe_mass_transfer(struct reb_simulation* const sim,
                                   struct rebx_operator*     const op,
                                   const double                     dt)
{
    struct rebx_extras* const rebx = sim->extras;

    /* indices of donor / accretor */
    int* donor_idx_ptr     = rebx_get_param(rebx, op->ap, "rlmt_donor");
    int* accretor_idx_ptr  = rebx_get_param(rebx, op->ap, "rlmt_accretor");
    if (!donor_idx_ptr || !accretor_idx_ptr){
        rebx_error(rebx,"Need to set rlmt_donor / rlmt_accretor.");  return;
    }
    int donor_idx    = *donor_idx_ptr;
    int accretor_idx = *accretor_idx_ptr;
    if (donor_idx>=sim->N || accretor_idx>=sim->N ||
        donor_idx<0     || accretor_idx<0){
        rebx_error(rebx,"rlmt_donor / rlmt_accretor out of range."); return;
    }

    struct reb_particle* d = &sim->particles[donor_idx];
    struct reb_particle* a = &sim->particles[accretor_idx];

    /* parameters ----------------------------------------------------------- */
    const double *Hp_ptr = rebx_get_param(rebx, d->ap, "rlmt_Hp");
    const double *mdot0_ptr = rebx_get_param(rebx, d->ap, "rlmt_mdot0");
    if (!Hp_ptr || !mdot0_ptr){
        rebx_error(rebx,"Set rlmt_Hp and rlmt_mdot0 on donor."); return;
    }
    const double *loss_ptr =
        rebx_get_param(rebx, op->ap, "rlmt_loss_fraction");
    const double *R_slope_ptr =
        rebx_get_param(rebx, d->ap, "rlmt_R_slope");
    const double *R_refM_ptr =
        rebx_get_param(rebx, d->ap, "rlmt_R_ref_mass");
    const double *R_refR_ptr =
        rebx_get_param(rebx, d->ap, "rlmt_R_ref_radius");
    /* CE parameters */
    const double *rho0_ptr  = rebx_get_param(rebx, op->ap, "ce_rho0");
    const double *arho_ptr  = rebx_get_param(rebx, op->ap, "ce_alpha_rho");
    const double *cs_ptr    = rebx_get_param(rebx, op->ap, "ce_cs");
    const double *acs_ptr   = rebx_get_param(rebx, op->ap, "ce_alpha_cs");
    const double *xmin_ptr  = rebx_get_param(rebx, op->ap, "ce_xmin");
    const double *Qd_ptr    = rebx_get_param(rebx, op->ap, "ce_Qd");
    const double *kick_ptr  = rebx_get_param(rebx, op->ap, "ce_max_kick_fraction");
    /* GW parameters */
    const double *c_ptr     = rebx_get_param(rebx, op->ap,"gw_c");
    const int    *gw_on_ptr = rebx_get_param(rebx, op->ap,"gw_decay_on");
    /* merge distance parameter (optional) */
    const double *merge_eps_ptr = rebx_get_param(rebx, op->ap, "merge_eps");

    /* convenient scalars */
    const double Hp     = *Hp_ptr;
    const double mdot0  = *mdot0_ptr;
    double loss_frac    =  loss_ptr ? *loss_ptr : 0.0;
    if (loss_frac<0) loss_frac=0; else if (loss_frac>1) loss_frac=1;

    /* ------------------------------------------------------------------ */
    /*  STEP 1  –  Roche‑lobe overflow (with AM correction, donor R(M))   */
    /* ------------------------------------------------------------------ */
    /* separation */
    const double dx = d->x - a->x;
    const double dy = d->y - a->y;
    const double dz = d->z - a->z;
    const double r  = sqrt(dx*dx + dy*dy + dz*dz);

    /* Eggleton Roche radius */
    const double q     = d->m / a->m;
    const double q13   = cbrt(q);
    const double RL    = r * (0.49*q13*q13) / (0.6*q13*q13 + log(1.+q13));

    /* Ritter exponential rate */
    const double Rd0   = d->r;
    double mdot        = -mdot0 * exp((Rd0-RL)/Hp);    /* <0 for loss */
    double dM          = mdot*dt;                      /* donor Δm (<0) */

    if (d->m + dM <= 0.)    dM = -d->m;               /* avoid neg mass */
    const double m_loss = -dM;                        /* positive */
    const double m_acc  = m_loss*(1.-loss_frac);

    /* masses after transfer */
    const double Md0 = d->m,  Ma0 = a->m;
    d->m -= m_loss;  a->m += m_acc;
    const double Md1 = d->m,  Ma1 = a->m;

    /* linear momentum: accreted gas arrives with donor velocity */
    const double px_a = Ma0*a->vx + m_acc*d->vx;
    const double py_a = Ma0*a->vy + m_acc*d->vy;
    const double pz_a = Ma0*a->vz + m_acc*d->vz;
    a->vx = px_a/Ma1;  a->vy=py_a/Ma1;  a->vz=pz_a/Ma1;

    /* angular‑momentum correction (same as previous version) */
    double Rcmx = (Md1*d->x + Ma1*a->x)/(Md1+Ma1);
    double Rcmy = (Md1*d->y + Ma1*a->y)/(Md1+Ma1);
    double Rcmz = (Md1*d->z + Ma1*a->z)/(Md1+Ma1);

    double rx_d=d->x-Rcmx, ry_d=d->y-Rcmy, rz_d=d->z-Rcmz;
    double rx_a=a->x-Rcmx, ry_a=a->y-Rcmy, rz_a=a->z-Rcmz;

    double Lcorr_x,Lcorr_y,Lcorr_z;
    cross3(rx_d,ry_d,rz_d, d->vx,d->vy,d->vz, &Lcorr_x,&Lcorr_y,&Lcorr_z);
    Lcorr_x*=m_acc; Lcorr_y*=m_acc; Lcorr_z*=m_acc;

    double Linj_x,Linj_y,Linj_z;
    cross3(rx_a,ry_a,rz_a, d->vx,d->vy,d->vz, &Linj_x,&Linj_y,&Linj_z);
    Linj_x*=m_acc;  Linj_y*=m_acc; Linj_z*=m_acc;

    double dLx=Lcorr_x-Linj_x, dLy=Lcorr_y-Linj_y, dLz=Lcorr_z-Linj_z;
    double r2_a = rx_a*rx_a+ry_a*ry_a+rz_a*rz_a;
    double dL2  = dLx*dLx+dLy*dLy+dLz*dLz;

    if (dL2>0 && r2_a>0){
        double tmpx,tmpy,tmpz;
        cross3(dLx,dLy,dLz, rx_a,ry_a,rz_a, &tmpx,&tmpy,&tmpz);
        double inv = 1./(Ma1*r2_a);
        double dvx=tmpx*inv, dvy=tmpy*inv, dvz=tmpz*inv;
        a->vx+=dvx; a->vy+=dvy; a->vz+=dvz;
        d->vx-=dvx*(Ma1/Md1); d->vy-=dvy*(Ma1/Md1); d->vz-=dvz*(Ma1/Md1);
    }

    /* donor R(M) update if desired */
    if (R_slope_ptr && *R_slope_ptr!=0.0){
        const double alpha=*R_slope_ptr;
        const double Mref = R_refM_ptr?*R_refM_ptr:Md0;
        const double Rref = R_refR_ptr?*R_refR_ptr:Rd0;
        d->r = Rref*pow(d->m/Mref,alpha);
    }
    const double Rd = d->r;

    /* ------------------------------------------------------------------ */
    /*  STEP 2 – Common‑envelope drag (unchanged except kick cap)         */
    /* ------------------------------------------------------------------ */
    if (rho0_ptr && cs_ptr && r<Rd){
        double rho0=*rho0_ptr, arho=arho_ptr?*arho_ptr:0.0;
        double cs0 =*cs_ptr , acs =acs_ptr ?*acs_ptr :0.0;
        double xmin=xmin_ptr?*xmin_ptr:1e-4;
        double Qd  =Qd_ptr?*Qd_ptr:0.0;
        double fmax=kick_ptr?*kick_ptr:0.3;

        double s=r/Rd;
        double rho=rho0*pow(s,arho);
        double cs =cs0 *pow(s,acs);

        double vrelx=a->vx-d->vx, vrely=a->vy-d->vy, vrelz=a->vz-d->vz;
        double vrel=sqrt(vrelx*vrelx+vrely*vrely+vrelz*vrelz);
        if (vrel>0 && rho>0){
            double I=I_prefactor(vrel/cs, xmin);
            double fc=4.*M_PI*sim->G*sim->G*a->m*rho/(vrel*vrel*vrel)*I;
            if (Qd>0 && a->r>0)
                fc+=M_PI*rho*a->r*a->r*vrel*Qd/a->m;
            double dvx=-fc*vrelx*dt, dvy=-fc*vrely*dt, dvz=-fc*vrelz*dt;
            double dv=sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
            double cap=fmax*vrel;
            if (dv>cap){
                double s=cap/dv; dvx*=s; dvy*=s; dvz*=s;
            }
            a->vx+=dvx; a->vy+=dvy; a->vz+=dvz;
            if (d->m>0){
                double f=a->m/d->m;
                d->vx-=dvx*f; d->vy-=dvy*f; d->vz-=dvz*f;
            }
        }
    }

    /* ------------------------------------------------------------------ */
    /*  STEP 3 – GW back‑reaction (unchanged)                             */
    /* ------------------------------------------------------------------ */
    if (c_ptr && gw_on_ptr && *gw_on_ptr && *c_ptr>0){
        struct reb_orbit o=reb_orbit_from_particle(sim->G,*d,*a);
        if (o.a>0){
            double m1=d->m, m2=a->m, c=*c_ptr, e=o.e, a0=o.a;
            double G3=sim->G*sim->G*sim->G;
            double da=-64./5.*G3*m1*m2*(m1+m2)/(pow(c,5)*pow(a0,3)*pow(1-e*e,3.5))
                      *(1+73./24.*e*e+37./96.*e*e*e*e);
            double de=-304./15.*e*G3*m1*m2*(m1+m2)/(pow(c,5)*pow(a0,4)*pow(1-e*e,2.5))
                      *(1+121./304.*e*e);
            double anew=a0+da*dt, enew=e+de*dt;
            if (anew<0)anew=0; if (enew<0)enew=0;
            if (enew>0.999999)enew=0.999999;

            double Mtot=m1+m2;
            double Rcmx=(m1*d->x+m2*a->x)/Mtot;
            double Rcmy=(m1*d->y+m2*a->y)/Mtot;
            double Rcmz=(m1*d->z+m2*a->z)/Mtot;
            double Vcmx=(m1*d->vx+m2*a->vx)/Mtot;
            double Vcmy=(m1*d->vy+m2*a->vy)/Mtot;
            double Vcmz=(m1*d->vz+m2*a->vz)/Mtot;

            struct reb_particle np=
                reb_particle_from_orbit(sim->G,*a,m1,
                                        anew,enew,
                                        o.inc,o.Omega,o.omega,o.f);

            double rx=np.x-a->x, ry=np.y-a->y, rz=np.z-a->z;
            double vx=np.vx-a->vx, vy=np.vy-a->vy, vz=np.vz-a->vz;

            double fd=m2/Mtot, fa=-m1/Mtot;
            d->x=Rcmx+fd*rx; d->y=Rcmy+fd*ry; d->z=Rcmz+fd*rz;
            a->x=Rcmx+fa*rx; a->y=Rcmy+fa*ry; a->z=Rcmz+fa*rz;
            d->vx=Vcmx+fd*vx; d->vy=Vcmy+fd*vy; d->vz=Vcmz+fd*vz;
            a->vx=Vcmx+fa*vx; a->vy=Vcmy+fa*vy; a->vz=Vcmz+fa*vz;
        }
    }

    /* ------------------------------------------------------------------ */
    /*  STEP 4 – merge if overlap                                          */
    /* ------------------------------------------------------------------ */
    double merge_eps = merge_eps_ptr ? *merge_eps_ptr : 1e-2; /* 1 cm */
    double dxm=d->x-a->x, dym=d->y-a->y, dzm=d->z-a->z;
    double r2merge=dxm*dxm+dym*dym+dzm*dzm;
    if (r2merge<=merge_eps*merge_eps){
        /* merge: keep donor index, remove accretor */
        double msum=d->m+a->m;
        if (msum<=0){                                              /* rare   */
            d->m=0;
        }else{
            /* barycentre position & velocity for linear‑momentum cons. */
            d->x=(d->m*d->x+a->m*a->x)/msum;
            d->y=(d->m*d->y+a->m*a->y)/msum;
            d->z=(d->m*d->z+a->m*a->z)/msum;
            d->vx=(d->m*d->vx+a->m*a->vx)/msum;
            d->vy=(d->m*d->vy+a->m*a->vy)/msum;
            d->vz=(d->m*d->vz+a->m*a->vz)/msum;
            d->m = msum;
            d->r = fmax(d->r, a->r);
        }
        /* remove accretor (higher index keeps donor pointer valid) */
        if (accretor_idx>donor_idx){
            reb_simulation_remove_particle(sim, accretor_idx, 1);
        }else{
            reb_simulation_remove_particle(sim, accretor_idx, 1);
            donor_idx--; /* donor index shifted */
            d=&sim->particles[donor_idx];
        }
        /* detach the operator if only one body left */
        if (sim->N<2){
            rebx_remove_operator(rebx, op);
            reb_simulation_move_to_com(sim);
            return;
        }
    }

    /* ------------------------------------------------------------------ */
    /*  STEP 5 – purge zero‑mass particles                                 */
    /* ------------------------------------------------------------------ */
    for (int i=sim->N-1;i>=0;--i){
        if (sim->particles[i].m<=0){
            reb_simulation_remove_particle(sim, i, 1);
            if (i==donor_idx || i==accretor_idx){
                /* if we removed one of the binary components, stop operator */
                rebx_remove_operator(rebx, op);
                reb_simulation_move_to_com(sim);
                return;
            }
        }
    }

    /* finally recentre */
    reb_simulation_move_to_com(sim);
}
