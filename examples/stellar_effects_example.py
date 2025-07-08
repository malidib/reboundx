import rebound
import reboundx
from reboundx import stellar_effects

sim = rebound.Simulation()
sim.add(m=1.)
sim.add(a=1.)
rebx = reboundx.Extras(sim)

stellar_effects.add_stellar_wind_mass_loss(rebx, body=0, mdot=1e-8)
stellar_effects.add_roche_lobe_mass_transfer(rebx, donor=0, accretor=1, mdot=1e-8)
stellar_effects.add_gravitational_wave_decay(rebx, body1=0, body2=1)

sim.integrate(10.)
print('masses:', sim.particles[0].m, sim.particles[1].m)
