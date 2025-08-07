import math
import unittest
import rebound
import reboundx

class TestStellarWinds(unittest.TestCase):
    def test_mass_loss_with_sse(self):
        sim = rebound.Simulation()
        sim.units = ('AU', 'yr', 'Msun')
        sim.G = 4 * math.pi**2
        sim.add(m=1.0, r=1.0)

        rebx = reboundx.Extras(sim)
        swml = rebx.load_operator('stellar_wind_mass_loss')
        rebx.add_operator(swml)
        tdw = rebx.load_operator('thermally_driven_winds')
        rebx.add_operator(tdw)
        sse = rebx.load_operator('stellar_evolution_sse')
        rebx.add_operator(sse)

        star = sim.particles[0]
        star.params['swml_eta'] = 0.5
        star.params['tdw_eta'] = 1.0
        star.params['tdw_T'] = 1.5e6

        initial_mass = star.m
        for _ in range(5):
            sim.integrate(sim.t + 2.0e4)

        self.assertLess(star.m, initial_mass)

if __name__ == '__main__':
    unittest.main()
