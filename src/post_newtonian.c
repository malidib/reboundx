/**
 * @file    post_newtonian.c
 * @brief   Post-Newtonian relativistic corrections up to 2.5 PN order.
 *
 * Implements the harmonic-coordinate point-mass equations of motion
 * from Kidder (1995, Phys. Rev. D 52, 821) for every massive pair.
 *  – 1 PN   : conservative peri-apsis precession
 *  – 1.5 PN : spin–orbit couplings
 *  – 2 PN   : conservative point-mass + spin–spin corrections
 *  – 2.5 PN : gravitational-wave radiation reaction
 *
 * Effect parameters
 * -----------------
 * c       (double, required)  – speed of light in simulation units
 * pn_1PN  (int,    optional)  – include 1 PN terms (default: 1)
 * pn_15PN (int,    optional)  – include 1.5 PN spin–orbit terms (default: 1)
 * pn_2PN  (int,    optional)  – include 2 PN terms (default: 1)
 * pn_25PN (int,    optional)  – include 2.5 PN terms (default: 1)
 *
 * Particle parameters
 * -------------------
 * pn_spin (reb_vec3d, optional) – spin angular momentum vector (mass·L²/T)
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"
typedef struct reb_vec3d reb_vec3d;

/* ------------------------------------------------------------------------- */
/* Pair-wise accelerations                                                   */
/* ------------------------------------------------------------------------- */
static inline void pn_add_pair(struct reb_simulation* const sim,
                               struct reb_particle*   const pi,
                               struct reb_particle*   const pj,
                               const double                 G,
                               const double                 c,
                               const int                    do1PN,
                               const int                    do15PN,
                               const int                    do2PN,
                               const int                    do25PN)
{
    /* Relative separation & velocity */
    const double dx  = pi->x  - pj->x;
    const double dy  = pi->y  - pj->y;
    const double dz  = pi->z  - pj->z;
    const double dvx = pi->vx - pj->vx;
    const double dvy = pi->vy - pj->vy;
    const double dvz = pi->vz - pj->vz;

    const double r2   = dx*dx + dy*dy + dz*dz;
    if (r2 == 0.0) return;
    const double r    = sqrt(r2);
    const double invr = 1.0/r;
    const double invr2= 1.0/r2;

    const double v2   = dvx*dvx + dvy*dvy + dvz*dvz;
    const double rdot = (dx*dvx + dy*dvy + dz*dvz)*invr;

    /* Mass factors */
    const double m   = pi->m + pj->m;
    const double mu  = pi->m * pj->m / m;
    const double eta = mu/m;

    /* Unit vector */
    const double nx = dx*invr;
    const double ny = dy*invr;
    const double nz = dz*invr;

    double ax = 0.0, ay = 0.0, az = 0.0;

    /* --- 1 PN conservative term (Kidder 1995 Eq.2.2b) --------------- */
    if (do1PN) {
        const double pref1 = G*m*invr2/(c*c);
        const double An = (1.0 + 3.0*eta)*v2
                        - 2.0*(2.0 + eta)*G*m*invr
                        - 1.5*eta*rdot*rdot;
        const double Bv = -2.0*(2.0 - eta)*rdot;
        ax += -pref1*(An*nx + Bv*dvx);
        ay += -pref1*(An*ny + Bv*dvy);
        az += -pref1*(An*nz + Bv*dvz);
    }

    /* --- 1.5 PN spin–orbit term (Kidder 1995 Eq.2.2c) ---------------- */
    if (do15PN) {
        struct rebx_extras* const rx = sim->extras;
        const struct reb_vec3d* Spi = rebx_get_param_vec(rx, pi->ap, "pn_spin");
        const struct reb_vec3d* Spj = rebx_get_param_vec(rx, pj->ap, "pn_spin");
        if (Spi || Spj) {
            /* retrieve spins, default zero */
            reb_vec3d Si = Spi ? *Spi : (reb_vec3d){0,0,0};
            reb_vec3d Sj = Spj ? *Spj : (reb_vec3d){0,0,0};

            /* combined & difference */
            double dm = pi->m - pj->m;
            double dm_m = dm/m;
            reb_vec3d S  = { Si.x + Sj.x, Si.y + Sj.y, Si.z + Sj.z };
            reb_vec3d D  = { m*(Sj.x/pj->m - Si.x/pi->m),
                             m*(Sj.y/pj->m - Si.y/pi->m),
                             m*(Sj.z/pj->m - Si.z/pi->m) };
            reb_vec3d A2 = { 2*S.x + dm_m*D.x,
                             2*S.y + dm_m*D.y,
                             2*S.z + dm_m*D.z };
            reb_vec3d A7 = { 7*S.x + 3*dm_m*D.x,
                             7*S.y + 3*dm_m*D.y,
                             7*S.z + 3*dm_m*D.z };
            reb_vec3d A3 = { 3*S.x + dm_m*D.x,
                             3*S.y + dm_m*D.y,
                             3*S.z + dm_m*D.z };

            /* helper cross products */
            reb_vec3d nxv = { ny*dvz - nz*dvy,
                              nz*dvx - nx*dvz,
                              nx*dvy - ny*dvx };
            double nxv_dot = nxv.x*A2.x + nxv.y*A2.y + nxv.z*A2.z;

            reb_vec3d vA7 = { dvy*A7.z - dvz*A7.y,
                              dvz*A7.x - dvx*A7.z,
                              dvx*A7.y - dvy*A7.x };

            reb_vec3d nA3 = { ny*A3.z - nz*A3.y,
                              nz*A3.x - nx*A3.z,
                              nx*A3.y - ny*A3.x };

            const double pref15 =  G*invr2*invr/(c*c*c);
            ax += pref15*( 6*nx*nxv_dot - vA7.x + 3*rdot*nA3.x );
            ay += pref15*( 6*ny*nxv_dot - vA7.y + 3*rdot*nA3.y );
            az += pref15*( 6*nz*nxv_dot - vA7.z + 3*rdot*nA3.z );
        }
    }

    /* --- 2 PN conservative + spin–spin (Kidder 1995 Eqs.2.2d,e) -------- */
    if (do2PN) {
        /* point-mass 2 PN */
        const double pref2 = G*m*invr2/(c*c*c*c);
        const double gmr   = G*m*invr;
        const double v4    = v2*v2;
        const double rdot2 = rdot*rdot;
        const double An2 = (3./4.)*(12.+29.*eta)*gmr*gmr
                         + eta*(3.-4.*eta)*v4
                         + (15./8.)*eta*(1.-3.*eta)*rdot2*rdot2
                         - (3./2.)*eta*(3.-4.*eta)*v2*rdot2
                         - (2.+25.*eta+2.*eta*eta)*gmr*rdot2;
        const double Bv2 = -0.5*rdot*( eta*(15.+4.*eta)*v2
                                    - (4.+41.*eta+8.*eta*eta)*gmr
                                    - 3.*eta*(3.+2.*eta)*rdot2 );
        ax += -pref2*(An2*nx + Bv2*dvx);
        ay += -pref2*(An2*ny + Bv2*dvy);
        az += -pref2*(An2*nz + Bv2*dvz);

        /* spin–spin coupling */
        struct rebx_extras* const rx = sim->extras;
        const struct reb_vec3d* Spi = rebx_get_param_vec(rx, pi->ap, "pn_spin");
        const struct reb_vec3d* Spj = rebx_get_param_vec(rx, pj->ap, "pn_spin");
        if (Spi && Spj) {
            reb_vec3d Si = *Spi, Sj = *Spj;
            double Sdot = Si.x*Sj.x + Si.y*Sj.y + Si.z*Sj.z;
            double nidSi = nx*Si.x + ny*Si.y + nz*Si.z;
            double nidSj = nx*Sj.x + ny*Sj.y + nz*Sj.z;
            const double prefSS = 3*G/(c*c*c*c)*invr2*invr2;
            double invm_i = 1.0/pi->m;
            double invm_j = 1.0/pj->m;
            double common = (Sdot - 5.*nidSi*nidSj)*invm_i*invm_j;
            ax += -prefSS*( common*nx + nidSj*Si.x*invm_i + nidSi*Sj.x*invm_j );
            ay += -prefSS*( common*ny + nidSj*Si.y*invm_i + nidSi*Sj.y*invm_j );
            az += -prefSS*( common*nz + nidSj*Si.z*invm_i + nidSi*Sj.z*invm_j );
        }
    }

    /* --- 2.5 PN radiation reaction (Kidder 1995 Eq.2.2f) -------------- */
    if (do25PN) {
        const double c2    = c*c;
        const double c4    = c2*c2;
        const double pref25 = 8.*G*G*m*m*eta/(5.*c4*c*r2*r);
        const double Cn = rdot*( 18.*v2 + (2./3.)*G*m*invr - 25.*rdot*rdot );
        const double Dv =        6.*v2 -  2.*G*m*invr - 15.*rdot*rdot;
        ax += pref25*( Cn*nx - Dv*dvx );
        ay += pref25*( Cn*ny - Dv*dvy );
        az += pref25*( Cn*nz - Dv*dvz );
    }

    /* --- Symmetric back-reaction on each body --------------------------- */
    const double fac_i = pj->m/m;
    const double fac_j = pi->m/m;
    pi->ax += fac_i*ax;  pj->ax -= fac_j*ax;
    pi->ay += fac_i*ay;  pj->ay -= fac_j*ay;
    pi->az += fac_i*az;  pj->az -= fac_j*az;
}

/* ------------------------------------------------------------------------- */
/* Force-kernel wrapper called by REBOUNDx                                   */
/* ------------------------------------------------------------------------- */
static void rebx_calculate_post_newtonian(struct reb_simulation* const sim,
                                          struct reb_particle*   const particles,
                                          const int                        N,
                                          const double                     c,
                                          const int                        do1PN,
                                          const int                        do15PN,
                                          const int                        do2PN,
                                          const int                        do25PN)
{
    const double G = sim->G;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            pn_add_pair(sim,
                        &particles[i], &particles[j],
                        G, c,
                        do1PN, do15PN, do2PN, do25PN);
        }
    }
}

void rebx_post_newtonian(struct reb_simulation* const sim,
                         struct rebx_force*      const force,
                         struct reb_particle*    const particles,
                         const int                            N)
{
    struct rebx_extras* const rx = sim->extras;
    const double* c_ptr = rebx_get_param(rx, force->ap, "c");
    if (!c_ptr || *c_ptr <= 0.0) {
        reb_simulation_error(sim,
            "post_newtonian: must supply a positive 'c'.");
        return;
    }

    /* default all PN orders on */
    int do1PN  = 1, do15PN = 1, do2PN = 1, do25PN = 1;
    int *ip;
    if ((ip = rebx_get_param(rx, force->ap, "pn_1PN" ))) do1PN  = *ip;
    if ((ip = rebx_get_param(rx, force->ap, "pn_15PN"))) do15PN = *ip;
    if ((ip = rebx_get_param(rx, force->ap, "pn_2PN" ))) do2PN  = *ip;
    if ((ip = rebx_get_param(rx, force->ap, "pn_25PN"))) do25PN = *ip;

    rebx_calculate_post_newtonian(sim, particles, N, *c_ptr,
                                  do1PN, do15PN, do2PN, do25PN);
}
