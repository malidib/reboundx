import rebound
import reboundx
import numpy as np

sim = rebound.Simulation()
sim.units = ('AU','yr','Msun')
#sim.G = 4*3.141592653589793**2
sim.integrator='whfast'
# donor star
sim.add(m=1.0, r=0.2)
# companion inside envelope
sim.add(m=0.81, a=0.19)
sim.particles[1].r = 0.001
sim.dt=sim.particles[1].P/100
sim.move_to_com()

rebx = reboundx.Extras(sim)
rl = rebx.load_operator('roche_lobe_mass_transfer')
rebx.add_operator(rl)

rl.params['rlmt_donor'] = 0
rl.params['rlmt_accretor'] = 1
sim.particles[0].params['rlmt_Hp'] = 0.01
sim.particles[0].params['rlmt_mdot0'] = 1e-5

# common-envelope drag parameters
rl.params['ce_rho0'] = 1e2
rl.params['ce_alpha_rho'] = -2.0
rl.params['ce_cs'] = 30.0
rl.params['ce_alpha_cs'] = 0.0
rl.params['ce_xmin'] = 0.01
rl.params['ce_Qd'] = 1.0

for t in np.logspace(-3,0,1000):
	sim.integrate(t)
	print(t,sim.particles[0].m, 'a={:.4f}'.format(sim.particles[1].a))
