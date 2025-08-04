/**
 * @file    post_newtonian.c
 * @brief   Post-Newtonian relativistic corrections up to 2.5PN order.
 *
 * This force adds the leading conservative 1PN precession terms and the
 * dissipative 2.5PN radiation–reaction terms for every pair of massive
 * bodies.  It can be combined with existing gravitational-wave decay
 * operators and other forces.
 *
 * The section after the dollar signs gets built into the documentation by
 * a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See
 * http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Post-Newtonian Relativistic Corrections$
 *
 * ======================= ===============================================
 * Authors                 REBOUNDx team
 * Based on                `Kidder 1995 <https://ui.adsabs.harvard.edu/abs/1995PhysRevD..52.821K/abstract>`_
 *                         `Peters 1964 <https://ui.adsabs.harvard.edu/abs/1964PhRv..136.1224P/abstract>`_
 * C Example               None
 * Python Example          None
 * ======================= ===============================================
 *
 * This force applies post-Newtonian accelerations through 2.5PN order for
 * each particle pair.  The 1PN term produces apsidal precession while the
 * 2.5PN term drives an inspiral due to gravitational-wave emission.
 *
 * **Effect Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * c (double)                   Yes         Speed of light in simulation units.
 * ============================ =========== ==================================================================
 *
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_post_newtonian(struct reb_simulation* const sim,
        struct reb_particle* const particles, const int N, const double c){
    const double G = sim->G;
    const double C2 = c*c;
    const double C5 = C2*c*c*c;

    for(int i=0;i<N;i++){
        for(int j=i+1;j<N;j++){
            struct reb_particle *pi = &particles[i];
            struct reb_particle *pj = &particles[j];

            const double dx = pi->x - pj->x;
            const double dy = pi->y - pj->y;
            const double dz = pi->z - pj->z;
            const double dvx = pi->vx - pj->vx;
            const double dvy = pi->vy - pj->vy;
            const double dvz = pi->vz - pj->vz;

            const double r2 = dx*dx + dy*dy + dz*dz;
            const double r = sqrt(r2);
            const double invr = 1./r;
            const double invr2 = 1./r2;
            const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
            const double rdot = (dx*dvx + dy*dvy + dz*dvz)*invr;

            const double m = pi->m + pj->m;
            const double eta = pi->m * pj->m / (m*m);

            const double nx = dx*invr;
            const double ny = dy*invr;
            const double nz = dz*invr;

            /* 1PN conservative term */
            const double A = (1.+3.*eta)*v2 - 2.*(2.+eta)*G*m*invr - 1.5*eta*rdot*rdot;
            const double B = 2.*(2.-eta)*rdot;
            const double pref1 = G*m*invr2/C2;
            double ax_rel = -pref1*(A*nx + B*dvx);
            double ay_rel = -pref1*(A*ny + B*dvy);
            double az_rel = -pref1*(A*nz + B*dvz);

            /* 2.5PN dissipative term */
            const double pref25 = 8.*G*G*m*m*eta/(5.*C5*r*r*r);
            const double Cterm = rdot*(3.*G*m*invr + v2);
            const double Dterm = G*m*invr + v2;
            ax_rel += pref25*(Cterm*nx - Dterm*dvx);
            ay_rel += pref25*(Cterm*ny - Dterm*dvy);
            az_rel += pref25*(Cterm*nz - Dterm*dvz);

            const double fac_i = pj->m/m;
            const double fac_j = pi->m/m;

            pi->ax += fac_i*ax_rel;
            pi->ay += fac_i*ay_rel;
            pi->az += fac_i*az_rel;
            pj->ax -= fac_j*ax_rel;
            pj->ay -= fac_j*ay_rel;
            pj->az -= fac_j*az_rel;
        }
    }
}

void rebx_post_newtonian(struct reb_simulation* const sim, struct rebx_force* const force,
        struct reb_particle* const particles, const int N){
    double* c_ptr = rebx_get_param(sim->extras, force->ap, "c");
    if(c_ptr == NULL){
        reb_simulation_error(sim, "Need to set speed of light in post_newtonian effect.\n");
        return;
    }
    rebx_calculate_post_newtonian(sim, particles, N, *c_ptr);
}

