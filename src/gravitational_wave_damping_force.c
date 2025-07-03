/*
 * @file    gravitational_wave_damping_force.c
 * @brief   Orbit-averaged gravitational-wave radiation reaction (2.5PN) force
 *          implemented through instantaneous accelerations.
 *
 * This applies gravitational-wave damping to binaries by modifying particle
 * accelerations such that the secular Peters (1964) evolution of the orbital
 * elements is reproduced.
 *
 * Effect Parameters
 * -----------------
 * gw_c (double)  [required]  Speed of light in the units used for the simulation.
 * coordinates (enum) [optional] Coordinate system for applying the force
 *                               (Jacobi, barycentric or particle). Defaults Jacobi.
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

/* Compute the GW acceleration on particle p relative to source. */
static struct reb_vec3d rebx_calculate_gw_force(struct reb_simulation* const sim,
                                               struct rebx_force* const force,
                                               struct reb_particle* p,
                                               struct reb_particle* source){
    struct reb_vec3d a = {0};
    struct rebx_extras* const rebx = sim->extras;

    const double* c_ptr = rebx_get_param(rebx, force->ap, "gw_c");
    if (c_ptr == NULL){
        rebx_error(rebx, "Need to set gw_c parameter for gravitational wave damping. See examples.\n");
        return a;
    }
    const double C = *c_ptr;
    const double G = sim->G;

    int err = 0;
    struct reb_orbit o = reb_orbit_from_particle_err(G, *p, *source, &err);
    if (err){
        return a;
    }
    if (o.e >= 1.0){
        return a;                       /* ignore hyperbolic orbits */
    }

    const double m1 = p->m;
    const double m2 = source->m;
    const double mu = m1*m2/(m1+m2);
    const double M  = m1 + m2;
    const double aorb = o.a;
    const double e    = o.e;

    const double fac = -64.0/5.0 * pow(G,3) * m1*m2*M / pow(C,5);
    const double e2 = e*e;
    const double denom = pow(aorb,3) * pow(1.0-e2,3.5);
    const double da_dt = fac/denom * (1.0 + 73.0/24.0*e2 + 37.0/96.0*e2*e2);
    const double de_dt = fac*(304.0/64.0)/ (pow(aorb,4)*pow(1.0-e2,2.5)) * e * (1.0 + 121.0/304.0*e2);

    const double dE_dt = G*m1*m2/(2.*aorb*aorb) * da_dt;
    const double prefL = sqrt(G*M*aorb*(1.0-e2));
    const double dL_dt = mu*(G*M/(2.*prefL))*((1.0-e2)*da_dt - 2.*aorb*e*de_dt);

    const double dx  = p->x - source->x;
    const double dy  = p->y - source->y;
    const double dz  = p->z - source->z;
    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;

    const double r  = sqrt(dx*dx + dy*dy + dz*dz);
    const double rdot = (dx*dvx + dy*dvy + dz*dvz)/r;
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
    const double vt2 = v2 - rdot*rdot;
    const double vt  = (vt2>0.) ? sqrt(vt2) : 0.;

    struct reb_vec3d nhat = { dx/r, dy/r, dz/r };
    struct reb_vec3d that = {0};
    if (vt > 0.){
        that.x = (dvx - rdot*nhat.x)/vt;
        that.y = (dvy - rdot*nhat.y)/vt;
        that.z = (dvz - rdot*nhat.z)/vt;
    }

    const double a_t = dL_dt/(mu*r);
    double a_r = 0.;
    if (fabs(rdot) > 0.){
        a_r = (dE_dt/mu - vt*a_t)/rdot;
    }

    a.x = a_r*nhat.x + a_t*that.x;
    a.y = a_r*nhat.y + a_t*that.y;
    a.z = a_r*nhat.z + a_t*that.z;

    return a;
}

void rebx_gravitational_wave_damping_force(struct reb_simulation* const sim,
                                           struct rebx_force* const force,
                                           struct reb_particle* const particles,
                                           const int N){
    int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI;
    if (ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive,
                   reference_name, rebx_calculate_gw_force, particles, N);
}

