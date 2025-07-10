import rebound
import reboundx

sim = rebound.Simulation()
sim.units = ('AU', 'yr', 'Msun')
#sim.G = 4*3.141592653589793**2
sim.ri_ias15.epsilon=0

sim.add(m=1.)  # accretor
sim.add(m=0.8, a=0.05)  # donor
sim.particles[1].r = 0.01

rebx = reboundx.Extras(sim)
rl = rebx.load_operator("roche_lobe_mass_transfer")
rebx.add_operator(rl)

rl.params['rlmt_donor'] = 1
rl.params['rlmt_accretor'] = 0
rl.params['rlmt_loss_fraction'] = 0.5
sim.particles[1].params['rlmt_Hp'] = 5e-2
sim.particles[1].params['rlmt_mdot0'] = 1e-2

sim.integrate(100.0)
print(f"m1={sim.particles[0].m}")
