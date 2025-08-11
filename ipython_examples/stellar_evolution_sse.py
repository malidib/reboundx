import rebound
import reboundx

sim = rebound.Simulation()
sim.units = ('AU','yr','Msun')
sim.G = 4*3.141592653589793**2
sim.add(m=1., r=1.)

rebx = reboundx.Extras(sim)
sse = rebx.load_operator("stellar_evolution_sse")
rebx.add_operator(sse)

print("Default main-sequence scaling:")
for i in range(5):
    sim.integrate(sim.t+1e3)
    R = sim.particles[0].r
    L = sim.particles[0].params['sse_L']
    print(f"t={sim.t:.0f} yr, R={R:.4f} Rsun, L={L:.4f} Lsun")

# Customize parameters

sim.particles[0].params['sse_R_coeff'] = 0.9
sim.particles[0].params['sse_R_exp'] = 0.6
sim.particles[0].params['sse_L_coeff'] = 1.2
sim.particles[0].params['sse_L_exp'] = 4.0

sim.t = 0.0


print("Default main-sequence evolution:")
for i in range(5):
    sim.integrate(sim.t+1e3)
    R = sim.particles[0].r
    L = sim.particles[0].params['sse_L']
    print(f"t={sim.t:.0f} yr, R={R:.4f} Rsun, L={L:.4f} Lsun")

# Customize parameters

sim.particles[0].params['sse_R_coeff'] = 0.9
sim.particles[0].params['sse_R_exp'] = 0.6
sim.particles[0].params['sse_L_coeff'] = 1.2
sim.particles[0].params['sse_L_exp'] = 4.0

sim.t = 0.0
sim.particles[0].params['sse_age'] = 0.0


print("\nCustom parameters:")
for i in range(5):
    sim.integrate(sim.t+1e3)
    R = sim.particles[0].r
    L = sim.particles[0].params['sse_L']
    print(f"t={sim.t:.0f} yr, R={R:.4f} Rsun, L={L:.4f} Lsun")
