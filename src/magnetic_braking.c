/**
 * @file    magnetic_braking.c
 * @brief   Spin down cool stars via magnetic braking.
 *
 * Implements a Verbunt-Zwaan / Kawaler type magnetic braking law
 * with saturation.  The braking torque is applied anti-parallel
 * to the current spin vector of any particle that has ``mb_on``
 * set in its additional parameters and a non-zero ``mb_convective``
 * flag.  The star's spin vector ``Omega`` is updated directly and
 * tidal coupling may subsequently transfer the removed angular
 * momentum to the orbit.
 *
 * **Particle parameters**
 *
 * ============================ =========== =====================================
 * Field (C type)               Required    Description
 * ============================ =========== =====================================
 * mb_on (int)                  Yes         Enable magnetic braking for particle.
 * mb_convective (int)          Yes         Flag indicating a convective envelope.
 * mb_omega_sat (double)        No          Saturation angular velocity.
 * I (double)                   Yes         Moment of inertia of the star.
 * Omega (reb_vec3d)            Yes         Spin angular frequency vector.
 * ============================ =========== =====================================
 *
 * **Operator parameters**
 *
 * ============================ =========== =====================================
 * Field (C type)               Required    Description
 * ============================ =========== =====================================
 * mb_K (double)                No          Braking constant K (default 2.7e47).
 * ============================ =========== =====================================
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_magnetic_braking(struct reb_simulation* const sim,
                            struct rebx_operator* const op,
                            const double dt){
    struct rebx_extras* const rebx = sim->extras;
    const int N = sim->N;  // operate over all particles

    const double* K_ptr = rebx_get_param(rebx, op->ap, "mb_K");
    const double K = K_ptr ? *K_ptr : 2.7e47;  // user-specified constant

    for(int i=0; i<N; i++){
        struct reb_particle* p = &sim->particles[i];

        const int* on_ptr  = rebx_get_param(rebx, p->ap, "mb_on");
        const int* conv_ptr = rebx_get_param(rebx, p->ap, "mb_convective");
        if(on_ptr == NULL || *on_ptr == 0) continue;
        if(conv_ptr == NULL || *conv_ptr == 0) continue;

        double* I_ptr = rebx_get_param(rebx, p->ap, "I");
        struct reb_vec3d* Omega = rebx_get_param(rebx, p->ap, "Omega");
        if(I_ptr == NULL || Omega == NULL) continue;
        if(*I_ptr <= 0.) continue;

        const double R = p->r;
        const double M = p->m;
        if(R <= 0. || M <= 0.) continue;

        const double omega = sqrt(Omega->x*Omega->x + Omega->y*Omega->y + Omega->z*Omega->z);
        if(!isfinite(omega) || omega == 0.) continue;

        const double* sat_ptr = rebx_get_param(rebx, p->ap, "mb_omega_sat");
        const double omega_sat = sat_ptr ? *sat_ptr : INFINITY;

        double omega_term = omega*omega*omega;
        if(omega > omega_sat){
            omega_term = omega * omega_sat * omega_sat;  // saturated regime
        }

        // Verbunt-Zwaan / Kawaler braking torque (anti-parallel to spin)
        const double torque = -K * pow(R,0.5) * pow(M,-0.5) * omega_term;

        const double fac = torque / (*I_ptr * omega); // dOmega/dt along -Omega_hat
        Omega->x += fac * Omega->x * dt;
        Omega->y += fac * Omega->y * dt;
        Omega->z += fac * Omega->z * dt;
    }
}

