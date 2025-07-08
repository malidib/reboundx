#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

/**
 * @file    gravitational_wave_decay.c
 * @brief   Simplified gravitational wave orbital decay (2.5PN) for binaries.
 *
 * $Relativity$
 *
 * ======================= ===============================================
 * Authors                 Codex
 * Implementation Paper    `Peters 1964 <https://ui.adsabs.harvard.edu/abs/1964PhRv..136.1224P/abstract>`_.
 * Based on                Peters (1964) orbit averaged decay rates.
 * C Example               :ref:`c_example_gw_orbital_decay`
 * Python Example          `GWOrbitalDecay.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GWOrbitalDecay.ipynb>`_.
 * ======================= ===============================================
 *
 * This operator applies orbital decay to particles orbiting the primary
 * due to gravitational wave emission using orbit-averaged formulae.
 * Only the first two particles are considered (primary and one companion).
 *
 * **Effect Parameters**
 *
 * ============================ =========== ========================================
 * Field (C type)               Required    Description
 * ============================ =========== ========================================
 * c (double)                   No          Speed of light (defaults to 299792458)
 * ============================ =========== ========================================
 */

static void gw_decay_update(struct reb_simulation* const sim, struct rebx_operator* const operator, double dt){
    if(sim->N<2){
        return; // need primary and companion
    }
    struct reb_particle* primary = &sim->particles[0];
    struct reb_particle* p = &sim->particles[1];
    double* cptr = rebx_get_param(sim->extras, operator->ap, "c");
    double c = cptr ? *cptr : 299792458.0; // m/s default
    double G = sim->G;
    int err = 0;
    struct reb_orbit o = reb_orbit_from_particle_err(sim->G, *p, *primary, &err);
    if(err){
        return;
    }
    double m1 = primary->m;
    double m2 = p->m;
    double a = o.a;
    double e = o.e;
    double coeff = (64./5.) * G*G*G * m1*m2*(m1+m2)/(c*c*c*c*c);
    double da_dt = -coeff/pow(a,3)/pow(1-e*e,3.5) * (1. + (73./24.)*e*e + (37./96.)*pow(e,4));
    double de_dt = -coeff*e/(pow(a,4)*pow(1-e*e,2.5)) * (304./15.) * (1. + (121./304.)*e*e);
    a += da_dt*dt;
    e += de_dt*dt;
    if(e<0) e=0;
    struct reb_particle newp = reb_particle_from_orbit(sim->G, *primary, m2, a, e, o.inc, o.Omega, o.omega, o.f);
    *p = newp;
}

void rebx_gw_orbital_decay(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    gw_decay_update(sim, operator, dt);
}
