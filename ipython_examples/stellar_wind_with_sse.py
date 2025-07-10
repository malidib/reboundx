import rebound
import reboundx

sim = rebound.Simulation()
sim.units = ('AU','yr','Msun')
sim.G = 4*3.141592653589793**2
sim.add(m=1., r=1.)

rebx = reboundx.Extras(sim)
wind = rebx.load_operator('stellar_wind_mass_loss')
rebx.add_operator(wind)

sse = rebx.load_operator('stellar_evolution_sse')
rebx.add_operator(sse)

# wind parameters
sim.particles[0].params['swml_eta'] = 0.5
sim.particles[0].params['swml_R'] = 1.0
sim.particles[0].params['swml_L'] = 1.0

print("Mass, radius and luminosity with winds:")
for i in range(5):
    sim.integrate(sim.t + 2e4)
    M = sim.particles[0].m
    R = sim.particles[0].r
    L = sim.particles[0].params['swml_L']
    print(f"t={sim.t:.0f} yr, M={M:.6f} Msun, R={R:.3f} Rsun, L={L:.3f} Lsun")
