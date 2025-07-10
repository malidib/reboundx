import rebound
import reboundx

sim = rebound.Simulation()
sim.units = ('AU','yr','Msun')
sim.G = 4*3.141592653589793**2
sim.add(m=1., r=1.)

rebx = reboundx.Extras(sim)
sse = rebx.load_operator("stellar_evolution_sse")
rebx.add_operator(sse)

sim.particles[0].params['sse_age'] = 0.0

for i in range(5):
    sim.integrate(sim.t+1e3)
    R = sim.particles[0].r
    L = sim.particles[0].params['swml_L']
    print(f"t={sim.t:.0f} yr, R={R:.4f} Rsun, L={L:.4f} Lsun")
