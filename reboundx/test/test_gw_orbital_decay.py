import rebound
import reboundx
import unittest

class TestGWDecay(unittest.TestCase):
    def test_load(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3, a=0.01, e=0.1)
        rebx = reboundx.Extras(sim)
        op = rebx.load_operator('gw_orbital_decay')
        rebx.add_operator(op)
        sim.dt = 1e-3
        sim.integrate(1.0)
        p1 = sim.particles[1]
        o = p1.orbit(sim.particles[0], sim.G)
        self.assertLess(o.a, 0.01)

if __name__ == '__main__':
    unittest.main()
