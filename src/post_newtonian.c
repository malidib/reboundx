/**
 * @file    post_newtonian.c
 * @brief   Post‑Newtonian relativistic corrections up to 2.5 PN order.
 *
 * Implements the harmonic‑coordinate point‑mass equations of motion
 * from Kidder (1995, Phys. Rev. D 52, 821) for every massive pair.
 *  – 1 PN   : conservative peri‑apsis precession
 *  – 2.5 PN : gravitational‑wave radiation reaction
 *
 * Effect parameters
 * -----------------
 * c  (double, required)  – speed of light in simulation units.
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"

/* ------------------------------------------------------------------------- */
/* Pair‑wise accelerations                                                   */
/* ------------------------------------------------------------------------- */
static inline void pn_add_pair(struct reb_particle* const pi,
                               struct reb_particle* const pj,
                               const double G, const double c)
{
    /* Relative separation and velocity */
    const double dx  = pi->x  - pj->x;
    const double dy  = pi->y  - pj->y;
    const double dz  = pi->z  - pj->z;
    const double dvx = pi->vx - pj->vx;
    const double dvy = pi->vy - pj->vy;
    const double dvz = pi->vz - pj->vz;

    const double r2   = dx*dx + dy*dy + dz*dz;
    if (r2 == 0.0) return;                     /* co‑located – ignore          */
    const double r    = sqrt(r2);
    const double invr = 1.0/r;
    const double invr2= 1.0/r2;

    const double v2   = dvx*dvx + dvy*dvy + dvz*dvz;
    const double rdot = (dx*dvx + dy*dvy + dz*dvz)*invr;

    /* Mass factors */
    const double m   = pi->m + pj->m;
    const double eta = (pi->m * pj->m)/(m*m);

    /* Unit vector n = r / |r| */
    const double nx = dx*invr;
    const double ny = dy*invr;
    const double nz = dz*invr;

    /* --- 1 PN conservative term (Kidder 1995, Eq. 2.2b) ------------------ */
    const double pref1 =  G*m*invr2/(c*c);                    /* G m / (r² c²)      */

    const double An =  (1.0 + 3.0*eta)*v2
                     - 2.0*(2.0 + eta)*G*m*invr
                     - 1.5*eta*rdot*rdot;

    const double Bv = -2.0*(2.0 - eta)*rdot;                  /* minus sign!        */

    double ax = -pref1*( An*nx + Bv*dvx );
    double ay = -pref1*( An*ny + Bv*dvy );
    double az = -pref1*( An*nz + Bv*dvz );

    /* --- 2.5 PN dissipative term (Kidder 1995, Eq. 2.2f) ----------------- */
    const double c2   = c*c;
    const double c4   = c2*c2;
    const double pref25 = 8.0*G*G*m*m*eta/(5.0*c4*c*r*r*r);   /* 8 G² m² η / (5 c⁵ r³) */

    const double Cn = rdot*( 18.0*v2 + (2.0/3.0)*G*m*invr - 25.0*rdot*rdot );
    const double Dv =        6.0*v2 -  2.0*G*m*invr - 15.0*rdot*rdot;

    ax += pref25*( Cn*nx - Dv*dvx );
    ay += pref25*( Cn*ny - Dv*dvy );
    az += pref25*( Cn*nz - Dv*dvz );

    /* --- Symmetric back‑reaction on each body --------------------------- */
    const double fac_i =  pj->m / m;   /* accel. experienced by i          */
    const double fac_j =  pi->m / m;

    pi->ax += fac_i*ax;   pj->ax -= fac_j*ax;
    pi->ay += fac_i*ay;   pj->ay -= fac_j*ay;
    pi->az += fac_i*az;   pj->az -= fac_j*az;
}

/* ------------------------------------------------------------------------- */
/* Force‑kernel wrapper called by REBOUNDx                                   */
/* ------------------------------------------------------------------------- */
static void rebx_calculate_post_newtonian(struct reb_simulation* const sim,
                                          struct reb_particle* const particles,
                                          const int N, const double c)
{
    const double G = sim->G;

    for (int i = 0; i < N; i++){
        for (int j = i + 1; j < N; j++){
            pn_add_pair(&particles[i], &particles[j], G, c);
        }
    }
}

void rebx_post_newtonian(struct reb_simulation* const sim,
                         struct rebx_force*    const force,
                         struct reb_particle*  const particles,
                         const int N)
{
    const double* c_ptr = rebx_get_param(sim->extras, force->ap, "c");
    if (c_ptr == NULL || *c_ptr <= 0.0){
        reb_simulation_error(sim,
            "post_newtonian: must supply a positive 'c' (speed of light).");
        return;
    }
    rebx_calculate_post_newtonian(sim, particles, N, *c_ptr);
}
