import rebound
import reboundx
import unittest
import numpy as np

class TestRocheLobeMassTransfer(unittest.TestCase):
    def test_mass_transfer_only(self):
        sim = rebound.Simulation()
        sim.units = ('AU', 'yr', 'Msun')
        sim.add(m=1.0)
        sim.add(m=0.8, a=0.05)
        sim.particles[1].r = 0.01
        rebx = reboundx.Extras(sim)
        rl = rebx.load_operator('roche_lobe_mass_transfer')
        rebx.add_operator(rl)
        rl.params['rlmt_donor'] = 1
        rl.params['rlmt_accretor'] = 0
        rl.params['rlmt_loss_fraction'] = 0.5
        sim.particles[1].params['rlmt_Hp'] = 5e-2
        sim.particles[1].params['rlmt_mdot0'] = 1e-2
        m0_a = sim.particles[0].m
        m0_d = sim.particles[1].m
        sim.integrate(20.0)
        self.assertGreater(sim.particles[0].m, m0_a)
        self.assertLess(sim.particles[1].m, m0_d)

    def test_common_envelope_drag(self):
        sim = rebound.Simulation()
        sim.units = ('AU', 'yr', 'Msun')
        sim.add(m=1.0, r=0.2)
        sim.add(m=0.01, a=0.1)
        sim.particles[1].r = 0.001
        sim.move_to_com()
        rebx = reboundx.Extras(sim)
        rl = rebx.load_operator('roche_lobe_mass_transfer')
        rebx.add_operator(rl)
        rl.params['rlmt_donor'] = 0
        rl.params['rlmt_accretor'] = 1
        sim.particles[0].params['rlmt_Hp'] = 0.01
        sim.particles[0].params['rlmt_mdot0'] = 0.0
        rl.params['ce_rho0'] = 1e-6
        rl.params['ce_alpha_rho'] = -2.0
        rl.params['ce_cs'] = 30.0
        rl.params['ce_alpha_cs'] = 0.0
        rl.params['ce_xmin'] = 0.01
        rl.params['ce_Qd'] = 1.0
        a0 = sim.particles[1].a
        sim.integrate(1.0)
        self.assertLess(sim.particles[1].a, a0)

    def test_gw_decay_only(self):
        sim = rebound.Simulation()
        sim.units = ('AU', 'yr', 'Msun')
        sim.add(m=1.4)
        sim.add(m=1.3, a=0.1, e=0.0)
        sim.particles[1].r = 0.001
        sim.ri_ias15.epsilon = 0
        rebx = reboundx.Extras(sim)
        rl = rebx.load_operator('roche_lobe_mass_transfer')
        rebx.add_operator(rl)
        rl.params['rlmt_donor'] = 1
        rl.params['rlmt_accretor'] = 0
        sim.particles[1].params['rlmt_Hp'] = 1e-5
        sim.particles[1].params['rlmt_mdot0'] = 0.0
        rl.params['gw_c'] = 63239.7263
        rl.params['gw_decay_on'] = 1
        sim.integrate(1.0)
        self.assertLess(sim.particles[1].a, 0.0)

    def test_mass_transfer_and_ce_merge(self):
        sim = rebound.Simulation()
        sim.units = ('AU', 'yr', 'Msun')
        sim.add(m=1.0, r=0.2)
        sim.add(m=0.81, a=0.19)
        sim.particles[1].r = 0.001
        sim.dt = sim.particles[1].P / 100
        sim.move_to_com()
        rebx = reboundx.Extras(sim)
        rl = rebx.load_operator('roche_lobe_mass_transfer')
        rebx.add_operator(rl)
        rl.params['rlmt_donor'] = 0
        rl.params['rlmt_accretor'] = 1
        sim.particles[0].params['rlmt_Hp'] = 0.01
        sim.particles[0].params['rlmt_mdot0'] = 1e-5
        rl.params['ce_rho0'] = 1e2
        rl.params['ce_alpha_rho'] = -2.0
        rl.params['ce_cs'] = 30.0
        rl.params['ce_alpha_cs'] = 0.0
        rl.params['ce_xmin'] = 0.01
        rl.params['ce_Qd'] = 1.0
        for t in np.logspace(-3, -1, 5):
            sim.integrate(t)
            if sim.N == 1:
                break
        self.assertEqual(sim.N, 1)

if __name__ == '__main__':
    unittest.main()
