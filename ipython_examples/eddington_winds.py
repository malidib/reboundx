import rebound
import reboundx

sim = rebound.Simulation()
sim.units = ('AU','yr','Msun')
sim.G = 4*3.141592653589793**2
sim.add(m=1.)

rebx = reboundx.Extras(sim)
wind = rebx.load_operator("eddington_winds")
rebx.add_operator(wind)

# Set luminosity above the Eddington limit
a = sim.particles[0]
a.params['eddw_L'] = 1e5

for i in range(500):
    sim.integrate(sim.t+2e4)
    print(f"t={sim.t:.1f} yr, M={a.m:.6f} Msun")
