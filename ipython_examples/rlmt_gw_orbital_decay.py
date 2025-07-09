import rebound
import reboundx
import numpy as np

sim = rebound.Simulation()
sim.units = ('AU','yr','Msun')
sim.G = 4*np.pi**2

sim.add(m=1.4)  # accretor
sim.add(m=1.3, a=0.01)  # donor
sim.particles[1].r = 0.003

rebx = reboundx.Extras(sim)
rl = rebx.load_operator("roche_lobe_mass_transfer")
rebx.add_operator(rl)

rl.params['rlmt_donor'] = 1
rl.params['rlmt_accretor'] = 0
sim.particles[1].params['rlmt_Hp'] = 1e-5
sim.particles[1].params['rlmt_mdot0'] = 1e-9
rl.params['gw_c'] = 63239.7263  # speed of light in AU/yr
rl.params["gw_decay_on"] = 1


for i in range(5):
    sim.integrate(sim.t+0.1)
    o = sim.particles[1].orbit(primary=sim.particles[0])
    print(f"t={sim.t:.2f} a={o.a:.6e} e={o.e:.6e}")

