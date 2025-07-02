/**
 * @file    gravitational_wave_damping.c
 * @brief   Orbit-averaged gravitational-wave radiation reaction (2.5PN) operator
 * @author  OpenAI Codex
 *
 * @section     LICENSE
 * Copyright (c) 2024 OpenAI Codex
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.
 * All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $General Relativity$
 *
 * ======================= ===============================================
 * Authors                 OpenAI Codex
 * Implementation Paper    `Peters 1964 <https://ui.adsabs.harvard.edu/abs/1964PhRv..136.1224P/abstract>`_
 * Based on                `Peters 1964 <https://ui.adsabs.harvard.edu/abs/1964PhRv..136.1224P/abstract>`_
 * ======================= ===============================================
 *
 * This implements the orbit-averaged 2.5PN gravitational-wave back reaction.
 * It damps the semimajor axis and eccentricity of binaries according to
 * the classical Peters (1964) formulae. The implementation updates osculating
 * orbital elements each timestep.
 *
 * **Effect Parameters**
 *
 * ============================ =========== ================================================
 * Field (C type)               Required    Description
 * ============================ =========== ================================================
 * gw_c (double)                Yes         Speed of light in the units used for the simulation.
 * coordinates (enum)           No          Type of elements to use for modification
 *                                         (Jacobi, barycentric or particle). Defaults Jacobi.
 * ============================ =========== ================================================
 *
 * **Particle Parameters**
 *
 * *None*
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

static struct reb_particle rebx_calculate_gw_damping(struct reb_simulation* const sim, struct rebx_operator* const operator, struct reb_particle* p, struct reb_particle* primary, const double dt){
    struct rebx_extras* const rebx = sim->extras;
    const double* const c_ptr = rebx_get_param(rebx, operator->ap, "gw_c");
    if (c_ptr == NULL){
        rebx_error(rebx, "Need to set gw_c parameter for gravitational wave damping. See examples.\n");
        return *p;
    }
    const double C = *c_ptr;
    int err = 0;
    struct reb_orbit o = reb_orbit_from_particle_err(sim->G, *p, *primary, &err);
    if(err){
        return *p;
    }
    const double m1 = p->m;
    const double m2 = primary->m;

    const double e = o.e;
    if (e >= 1.0){
        return *p; // hyperbolic orbits ignored
    }

    const double G = sim->G;
    const double a = o.a;

    const double fac = -64.0/5.0 * pow(G,3) * m1*m2*(m1+m2) / (pow(C,5));
    const double one_minus_e2 = 1.0 - e*e;
    const double e2 = e*e;
    const double da_dt = fac / (pow(a,3) * pow(one_minus_e2,3.5)) * (1.0 + 73.0/24.0*e2 + 37.0/96.0*e2*e2);
    const double de_dt = fac * (304.0/64.0) / (pow(a,4) * pow(one_minus_e2,2.5)) * e * (1.0 + 121.0/304.0*e2);

    double anew = a + da_dt*dt;
    double enew = e + de_dt*dt;

    if(anew < 0.){ // avoid negative semimajor axis
        anew = 0.;
    }
    if(enew < 0.) enew = 0.;
    if(enew > 0.999999) enew = 0.999999;

    struct reb_particle newp = reb_particle_from_orbit(G, *primary, m1, anew, enew, o.inc, o.Omega, o.omega, o.f);
    return newp;
}

void rebx_gravitational_wave_damping(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    int* ptr = rebx_get_param(sim->extras, operator->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if(ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_tools_com_ptm(sim, operator, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_gw_damping, dt);
}

