/**
 * @file    post_newtonian.c
 * @brief   Post-Newtonian relativistic corrections: 2PN (PM+SS), 2.5PN (RR).
 *
 * Implements the harmonic-coordinate point-mass equations of motion from
 * Kidder (1995, Phys. Rev. D 52, 821) for every massive pair, *omitting 1PN*.
 *
 *  – 2 PN   : conservative point-mass (PM) + spin–spin (SS) corrections
 *  – 2.5 PN : gravitational-wave radiation reaction (RR)
 *
 * Effect parameters
 * -----------------
 * c        (double, required)  – speed of light in simulation units
 * pn_2PN   (double,bool, opt.) – include 2 PN terms (PM + SS)   (default: 1)
 * pn_25PN  (double,bool, opt.) – include 2.5 PN terms           (default: 1)
 *
 * Particle parameters
 * -------------------
 * pn_spin (reb_vec3d, optional) – spin angular momentum vector (mass·L^2/T)
 *
 * Notes
 * -----
 *  • Spins must be supplied as *physical* angular momenta. For dimensionless
 *    spins χ_i, set S_i = χ_i * (G m_i^2 / c) in the simulation's unit system.
 *  • This effect updates *accelerations only*; spin precession ODEs are not
 *    included here.
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"

typedef struct reb_vec3d reb_vec3d;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------------- */
/* Pair-wise relative acceleration builder                                   */
/* ------------------------------------------------------------------------- */
static inline void pn_add_pair(struct reb_simulation* const sim,
                               struct reb_particle*   const pi,
                               struct reb_particle*   const pj,
                               const double                 G,
                               const double                 c,
                               const int                    do2PN,
                               const int                    do25PN)
{
    /* Mass guards: PN is defined for massive pairs only */
    if (!(pi->m > 0.0 && pj->m > 0.0)) return;

    /* Relative separation & velocity */
    const double dx  = pi->x  - pj->x;
    const double dy  = pi->y  - pj->y;
    const double dz  = pi->z  - pj->z;
    const double dvx = pi->vx - pj->vx;
    const double dvy = pi->vy - pj->vy;
    const double dvz = pi->vz - pj->vz;

    const double r2 = dx*dx + dy*dy + dz*dz;
    if (!(r2 > 0.0)) return;                           /* coincident or NaN */
    const double r     = sqrt(r2);
    const double invr  = 1.0/r;
    const double invr2 = 1.0/r2;

    const double v2   = dvx*dvx + dvy*dvy + dvz*dvz;
    const double rdot = (dx*dvx + dy*dvy + dz*dvz) * invr;

    /* Mass combinations */
    const double m   = pi->m + pj->m;
    const double mu  = (pi->m * pj->m) / m;           /* reduced mass */
    const double eta = mu / m;

    /* Unit separation vector n */
    const double nx = dx*invr;
    const double ny = dy*invr;
    const double nz = dz*invr;

    /* Build relative-acceleration contributions */
    double ax = 0.0, ay = 0.0, az = 0.0;

    /* --------------------------------------------------------------------- */
    /* 2 PN conservative: point‑mass (Kidder 1995 Eq. 2.2d)                  */
    /* --------------------------------------------------------------------- */
    if (do2PN){
        const double gmr    = G*m*invr;                 /* (speed)^2 */
        const double rdot2  = rdot*rdot;
        const double v4     = v2*v2;
        const double pref2  = (G*m*invr2)/(c*c*c*c);    /* (G m / r^2) / c^4 */

        double A2 = (3.0/4.0)*(12.0 + 29.0*eta)*gmr*gmr
                  + eta*(3.0 - 4.0*eta)*v4
                  + (15.0/8.0)*eta*(1.0 - 3.0*eta)*rdot2*rdot2
                  - (3.0/2.0)*eta*(3.0 - 4.0*eta)*v2*rdot2
                  - (2.0 + 25.0*eta + 2.0*eta*eta)*gmr*rdot2
                  - 0.5*eta*(13.0 - 4.0*eta)*gmr*v2;   /* missing term restored */

        const double Bv2 = -0.5*rdot*( eta*(15.0 + 4.0*eta)*v2
                                     - (4.0 + 41.0*eta + 8.0*eta*eta)*gmr
                                     - 3.0*eta*(3.0 + 2.0*eta)*rdot2 );

        ax += -pref2 * ( A2*nx + Bv2*dvx );
        ay += -pref2 * ( A2*ny + Bv2*dvy );
        az += -pref2 * ( A2*nz + Bv2*dvz );
    }

    /* --------------------------------------------------------------------- */
    /* 2 PN spin–spin (Kidder 1995 Eq. 2.2e), physical spins.                */
    /* Prefactor: 3 G/(μ c^2 r^4) times the standard SS bracket.             */
    /* --------------------------------------------------------------------- */
    if (do2PN){
        struct rebx_extras* const rx = sim->extras;
        const reb_vec3d* Spi = rebx_get_param_vec(rx, pi->ap, "pn_spin");
        const reb_vec3d* Spj = rebx_get_param_vec(rx, pj->ap, "pn_spin");
        if (Spi && Spj){
            const reb_vec3d Si = *Spi;
            const reb_vec3d Sj = *Spj;

            const double Sdot   = Si.x*Sj.x + Si.y*Sj.y + Si.z*Sj.z;
            const double nidSi  = nx*Si.x + ny*Si.y + nz*Si.z;
            const double nidSj  = nx*Sj.x + ny*Sj.y + nz*Sj.z;

            const double prefSS = (3.0*G)/(mu * c*c) * invr2 * invr2;   /* 3G/(μ c^2 r^4) */

            ax += -prefSS * ( (Sdot - 5.0*nidSi*nidSj)*nx + nidSj*Si.x + nidSi*Sj.x );
            ay += -prefSS * ( (Sdot - 5.0*nidSi*nidSj)*ny + nidSj*Si.y + nidSi*Sj.y );
            az += -prefSS * ( (Sdot - 5.0*nidSi*nidSj)*nz + nidSj*Si.z + nidSi*Sj.z );
        }
    }

    /* --------------------------------------------------------------------- */
    /* 2.5 PN radiation reaction (Kidder 1995 Eq. 2.2f)                      */
    /* --------------------------------------------------------------------- */
    if (do25PN){
        const double c2 = c*c;
        const double c4 = c2*c2;
        const double pref25 = (8.0/(5.0)) * (G*G) * m*m * eta / (c4*c) * invr2 * invr;   /* 8 G^2 m^2 η /(5 c^5 r^3) */

        const double Cn = rdot * ( 18.0*v2 + (2.0/3.0)*(G*m*invr) - 25.0*rdot*rdot );
        const double Dv =          6.0*v2 -  2.0*(G*m*invr)      - 15.0*rdot*rdot;

        ax += pref25 * ( Cn*nx - Dv*dvx );
        ay += pref25 * ( Cn*ny - Dv*dvy );
        az += pref25 * ( Cn*nz - Dv*dvz );
    }

    /* --------------------------------------------------------------------- */
    /* Symmetric back‑reaction on each body (maps relative accel to bodies)   */
    /* --------------------------------------------------------------------- */
    const double fac_i = pj->m / m;
    const double fac_j = pi->m / m;

    pi->ax +=  fac_i * ax;  pj->ax -=  fac_j * ax;
    pi->ay +=  fac_i * ay;  pj->ay -=  fac_j * ay;
    pi->az +=  fac_i * az;  pj->az -=  fac_j * az;
}

/* ------------------------------------------------------------------------- */
/* Force-kernel wrapper called by REBOUNDx                                   */
/* ------------------------------------------------------------------------- */
static void rebx_calculate_post_newtonian(struct reb_simulation* const sim,
                                          struct reb_particle*   const particles,
                                          const int                        N,
                                          const double                     c,
                                          const int                        do2PN,
                                          const int                        do25PN)
{
    const double G = sim->G;
    if (!(isfinite(G) && isfinite(c) && c > 0.0)) return;

    for (int i = 0; i < N; i++){
        for (int j = i+1; j < N; j++){
            pn_add_pair(sim, &particles[i], &particles[j],
                        G, c, do2PN, do25PN);
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Public entry point for REBOUNDx                                           */
/* ------------------------------------------------------------------------- */
void rebx_post_newtonian(struct reb_simulation* const sim,
                         struct rebx_force*      const force,
                         struct reb_particle*    const particles,
                         const int                            N)
{
    struct rebx_extras* const rx = sim->extras;

    const double* c_ptr = rebx_get_param(rx, force->ap, "c");
    if (!(c_ptr && isfinite(*c_ptr) && *c_ptr > 0.0)){
        reb_simulation_error(sim, "post_newtonian: must supply a positive 'c' (speed of light).");
        return;
    }

    int do2PN = 1, do25PN = 1;

    const double* d;

    d = rebx_get_param(rx, force->ap, "pn_2PN");
    if (d) do2PN = (*d != 0.0);

    d = rebx_get_param(rx, force->ap, "pn_25PN");
    if (d) do25PN = (*d != 0.0);

    rebx_calculate_post_newtonian(sim, particles, N, *c_ptr, do2PN, do25PN);
}
