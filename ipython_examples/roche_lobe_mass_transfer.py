import rebound
import reboundx

sim = rebound.Simulation()
sim.units = ('AU', 'yr', 'Msun')
sim.G = 4*3.141592653589793**2

sim.add(m=1.)  # accretor
sim.add(m=0.8, a=0.05)  # donor
sim.particles[1].r = 0.01

rebx = reboundx.Extras(sim)
rl = rebx.load_operator("roche_lobe_mass_transfer")
rebx.add_operator(rl)

rl.params['rlmt_donor'] = 1
rl.params['rlmt_accretor'] = 0
sim.particles[1].params['rlmt_Hp'] = 5e-4
sim.particles[1].params['rlmt_mdot0'] = 1e-5

sim.integrate(1.0)
print(f"m1={sim.particles[0].m} m2={sim.particles[1].m}")
