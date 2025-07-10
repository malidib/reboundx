import rebound
import reboundx

sim = rebound.Simulation()
sim.units = ('AU','yr','Msun')
sim.G = 4*3.141592653589793**2
sim.add(m=1.)

rebx = reboundx.Extras(sim)
wind = rebx.load_operator("stellar_wind_mass_loss")
rebx.add_operator(wind)

sim.particles[0].params['swml_eta'] = 0.5
sim.particles[0].params['swml_R'] = 1.0
sim.particles[0].params['swml_L'] = 1.0

for i in range(5):
    sim.integrate(sim.t+2e4)
    print(f"t={sim.t:.1f} yr, M={sim.particles[0].m:.6f} Msun")
