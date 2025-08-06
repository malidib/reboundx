/**
 * @file    post_newtonian.c
 * @brief   Post‑Newtonian relativistic corrections up to 2.5 PN order.
 *
 * Implements the harmonic‑coordinate point‑mass equations of motion
 * from Kidder (1995, Phys. Rev. D 52, 821) for every massive pair.
 *  – 1 PN   : conservative peri‑apsis precession
 *  – 1.5 PN : spin–orbit couplings
 *  – 2 PN   : higher‑order conservative and spin–spin corrections
 *  – 2.5 PN : gravitational‑wave radiation reaction
 *
 * Effect parameters
 * -----------------
 * c         (double, required) – speed of light in simulation units.
 * pn_1PN    (int,    optional) – include 1 PN terms (default 1).
 * pn_15PN   (int,    optional) – include 1.5 PN spin–orbit terms (default 1).
 * pn_2PN    (int,    optional) – include 2 PN conservative/spin–spin terms (default 1).
 * pn_25PN   (int,    optional) – include 2.5 PN radiation‑reaction terms (default 1).
 *
 * Particle parameters
 * -------------------
 * pn_spin   (reb_vec3d, optional) – spin angular momentum vector of the particle
 *                                    in simulation units (mass·length²/time).
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"

/* ------------------------------------------------------------------------- */
/* Pair‑wise accelerations                                                   */
/* ------------------------------------------------------------------------- */
static inline void pn_add_pair(struct reb_simulation* const sim,
                               struct reb_particle* const pi,
                               struct reb_particle* const pj,
                               const double G, const double c,
                               const int do1PN, const int do15PN,
                               const int do2PN, const int do25PN)
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

    double ax = 0.0, ay = 0.0, az = 0.0;

    /* --- 1 PN conservative term (Kidder 1995, Eq. 2.2b) ------------------ */
    if (do1PN){
        const double pref1 =  G*m*invr2/(c*c);                /* G m / (r² c²)      */

        const double An =  (1.0 + 3.0*eta)*v2
                         - 2.0*(2.0 + eta)*G*m*invr
                         - 1.5*eta*rdot*rdot;

        const double Bv = -2.0*(2.0 - eta)*rdot;              /* minus sign!        */

        ax += -pref1*( An*nx + Bv*dvx );
        ay += -pref1*( An*ny + Bv*dvy );
        az += -pref1*( An*nz + Bv*dvz );
    }

    /* --- 1.5 PN spin–orbit term (Kidder 1995, Eq. 2.2c) ----------------- */
    if (do15PN){
        const struct reb_vec3d* Si_ptr = rebx_get_param_vec(sim->extras, pi->ap, "pn_spin");
        const struct reb_vec3d* Sj_ptr = rebx_get_param_vec(sim->extras, pj->ap, "pn_spin");
        struct reb_vec3d Si = {0.,0.,0.};
        struct reb_vec3d Sj = {0.,0.,0.};
        if (Si_ptr) Si = *Si_ptr;
        if (Sj_ptr) Sj = *Sj_ptr;

        const double dm = pi->m - pj->m;
        const double dm_over_m = dm/m;

        const struct reb_vec3d S = { Si.x + Sj.x, Si.y + Sj.y, Si.z + Sj.z };
        const struct reb_vec3d Delta = {
            m*(Sj.x/pj->m - Si.x/pi->m),
            m*(Sj.y/pj->m - Si.y/pi->m),
            m*(Sj.z/pj->m - Si.z/pi->m)
        };

        const struct reb_vec3d term2S = {
            2.*S.x + dm_over_m*Delta.x,
            2.*S.y + dm_over_m*Delta.y,
            2.*S.z + dm_over_m*Delta.z
        };
        const struct reb_vec3d term7S = {
            7.*S.x + 3.*dm_over_m*Delta.x,
            7.*S.y + 3.*dm_over_m*Delta.y,
            7.*S.z + 3.*dm_over_m*Delta.z
        };
        const struct reb_vec3d term3S = {
            3.*S.x + dm_over_m*Delta.x,
            3.*S.y + dm_over_m*Delta.y,
            3.*S.z + dm_over_m*Delta.z
        };

        /* cross products */
        const struct reb_vec3d nxv = { ny*dvz - nz*dvy,
                                       nz*dvx - nx*dvz,
                                       nx*dvy - ny*dvx };
        const double nxv_dot = nxv.x*term2S.x + nxv.y*term2S.y + nxv.z*term2S.z;

        const struct reb_vec3d vxt = { dvy*term7S.z - dvz*term7S.y,
                                       dvz*term7S.x - dvx*term7S.z,
                                       dvx*term7S.y - dvy*term7S.x };

        const struct reb_vec3d nxs = { ny*term3S.z - nz*term3S.y,
                                       nz*term3S.x - nx*term3S.z,
                                       nx*term3S.y - ny*term3S.x };

        const double pref15 = G*invr2*invr/(c*c*c);           /* G / (c³ r³)        */

        ax += pref15*( 6.*nx*nxv_dot - vxt.x + 3.*rdot*nxs.x );
        ay += pref15*( 6.*ny*nxv_dot - vxt.y + 3.*rdot*nxs.y );
        az += pref15*( 6.*nz*nxv_dot - vxt.z + 3.*rdot*nxs.z );
    }

    /* --- 2 PN conservative and spin–spin terms (Kidder 1995, Eqs. 2.2d–e) */
    if (do2PN){
        const double c2 = c*c;
        const double pref2 = G*m*invr2/(c2*c2);               /* G m / (r² c⁴)      */

        const double gmr = G*m*invr;                          /* G m / r            */

        const double An2 = (3./4.)*(12.+29.*eta)*gmr*gmr
                         + eta*(3.-4.*eta)*v2*v2
                         + (15./8.)*eta*(1.-3.*eta)*rdot*rdot*rdot*rdot
                         - (3./2.)*eta*(3.-4.*eta)*v2*rdot*rdot
                         - 0.5*eta*(13.-4.*eta)*gmr*v2
                         - (2.+25.*eta+2.*eta*eta)*gmr*rdot*rdot;

        const double Bv2 = -0.5*rdot*( eta*(15.+4.*eta)*v2
                                     - (4.+41.*eta+8.*eta*eta)*gmr
                                     - 3.*eta*(3.+2.*eta)*rdot*rdot );

        ax += -pref2*( An2*nx + Bv2*dvx );
        ay += -pref2*( An2*ny + Bv2*dvy );
        az += -pref2*( An2*nz + Bv2*dvz );

        const struct reb_vec3d* Si_ptr = rebx_get_param_vec(sim->extras, pi->ap, "pn_spin");
        const struct reb_vec3d* Sj_ptr = rebx_get_param_vec(sim->extras, pj->ap, "pn_spin");
        struct reb_vec3d Si = {0.,0.,0.};
        struct reb_vec3d Sj = {0.,0.,0.};
        if (Si_ptr) Si = *Si_ptr;
        if (Sj_ptr) Sj = *Sj_ptr;

        const double mu = pi->m*pj->m/m;
        const double S1dotS2 = Si.x*Sj.x + Si.y*Sj.y + Si.z*Sj.z;
        const double ndotS1 = nx*Si.x + ny*Si.y + nz*Si.z;
        const double ndotS2 = nx*Sj.x + ny*Sj.y + nz*Sj.z;

        const double prefSS = 3.*G*G/(c2*c2*mu)*invr2*invr2;

        const double ssx = nx*S1dotS2 + Si.x*ndotS2 + Sj.x*ndotS1 - 5.*nx*ndotS1*ndotS2;
        const double ssy = ny*S1dotS2 + Si.y*ndotS2 + Sj.y*ndotS1 - 5.*ny*ndotS1*ndotS2;
        const double ssz = nz*S1dotS2 + Si.z*ndotS2 + Sj.z*ndotS1 - 5.*nz*ndotS1*ndotS2;

        ax += -prefSS*ssx;
        ay += -prefSS*ssy;
        az += -prefSS*ssz;
    }

    /* --- 2.5 PN dissipative term (Kidder 1995, Eq. 2.2f) ----------------- */
    if (do25PN){
        const double c2   = c*c;
        const double c4   = c2*c2;
        const double pref25 = 8.0*G*G*m*m*eta/(5.0*c4*c*r*r*r); /* 8 G² m² η / (5 c⁵ r³) */

        const double Cn = rdot*( 18.0*v2 + (2.0/3.0)*G*m*invr - 25.0*rdot*rdot );
        const double Dv =        6.0*v2 -  2.0*G*m*invr - 15.0*rdot*rdot;

        ax += pref25*( Cn*nx - Dv*dvx );
        ay += pref25*( Cn*ny - Dv*dvy );
        az += pref25*( Cn*nz - Dv*dvz );
    }

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
                                          const int N, const double c,
                                          const int do1PN, const int do15PN,
                                          const int do2PN, const int do25PN)
{
    const double G = sim->G;

    for (int i = 0; i < N; i++){
        for (int j = i + 1; j < N; j++){
            pn_add_pair(sim, &particles[i], &particles[j], G, c,
                        do1PN, do15PN, do2PN, do25PN);
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
    int do1PN  = 1;
    int do15PN = 1;
    int do2PN  = 1;
    int do25PN = 1;
    int* ptr;
    ptr = rebx_get_param(sim->extras, force->ap, "pn_1PN");
    if (ptr) do1PN = *ptr;
    ptr = rebx_get_param(sim->extras, force->ap, "pn_15PN");
    if (ptr) do15PN = *ptr;
    ptr = rebx_get_param(sim->extras, force->ap, "pn_2PN");
    if (ptr) do2PN = *ptr;
    ptr = rebx_get_param(sim->extras, force->ap, "pn_25PN");
    if (ptr) do25PN = *ptr;
    rebx_calculate_post_newtonian(sim, particles, N, *c_ptr,
                                  do1PN, do15PN, do2PN, do25PN);
}
