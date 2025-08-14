import rebound
import reboundx

def test_eddington_mass_loss():
    sim = rebound.Simulation()
    sim.units = ('AU','yr','Msun')
    sim.G = 4*3.141592653589793**2
    sim.add(m=1.)

    rebx = reboundx.Extras(sim)
    op = rebx.load_operator('eddington_winds')
    rebx.add_operator(op)

    star = sim.particles[0]
    star.params['sse_L'] = 1e5
    m0 = star.m
    sim.integrate(1e3)
    assert star.m < m0
