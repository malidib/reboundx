#include "rebound.h"
#include "reboundx.h"
#include <stdio.h>

int main(){
    struct reb_simulation* sim = reb_create_simulation();
    reb_simulation_init(sim);
    reb_add_fmt(sim, "m", 1.0); // primary
    reb_add_fmt(sim, "m a e", 1e-3, 0.01, 0.1); // companion
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_operator* gw = rebx_load_operator(rebx, "gw_orbital_decay");
    rebx_add_operator(rebx, gw);
    sim->dt = 1e-3;
    for(int i=0;i<1000;i++){
        reb_integrate(sim, sim->t + sim->dt);
    }
    printf("final a = %e\n", reb_orbit(sim, sim->particles[1], sim->particles[0]).a);
    rebx_free(rebx);
    reb_free_simulation(sim);
}
