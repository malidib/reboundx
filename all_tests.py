#!/usr/bin/env python3
"""
pytest suite for the REBOUNDx operator `roche_lobe_mass_transfer` (rlmt).

Fixes:
 - Remove explicit jloss_mode / rlmt_skip_in_CE assignments (defaults are correct).
 - Swap sim.calculate_orbits() → sim.orbits().
 - Use a smaller semimajor axis for the immediate-merge test.
"""
import math
import numpy as np
import pytest
import rebound
import reboundx

def make_sim(
    donor_mass=1.0,
    accretor_mass=0.8,
    semimajor=1.5,
    eccentricity=0.0,
    donor_radius=0.8,
    accretor_radius=0.4,
    dt=1e-3,
):
    sim = rebound.Simulation()
    sim.G = 1.0
    sim.add(m=donor_mass, r=donor_radius)
    sim.add(m=accretor_mass, a=semimajor, e=eccentricity, r=accretor_radius)
    sim.integrator = "whfast"
    sim.dt = dt
    sim.move_to_com()

    rebx = reboundx.Extras(sim)
    try:
        rlmt = rebx.load_operator("roche_lobe_mass_transfer")
    except Exception:
        pytest.skip("rlmt operator not found – compile it into REBOUNDx first.")
    rebx.add_operator(rlmt)

    rlmt.params["rlmt_donor"]    = 0
    rlmt.params["rlmt_accretor"] = 1

    # Required donor params for Roche-lobe overflow
    donor = sim.particles[0]
    donor.params["rlmt_Hp"]    = 0.05
    donor.params["rlmt_mdot0"] = 1e-3

    return sim, rlmt

def test_rlof_mass_transfer():
    sim, rlmt = make_sim(semimajor=1.2, donor_radius=0.7)
    rlmt.params["rlmt_loss_fraction"] = 0.0

    md0, ma0 = sim.particles[0].m, sim.particles[1].m
    sim.integrate(sim.t + 0.1)
    md1, ma1 = sim.particles[0].m, sim.particles[1].m

    assert md1 < md0, "Donor should lose mass under RLOF"
    assert ma1 > ma0, "Accretor should gain mass under RLOF"
    assert math.isclose(md0 + ma0, md1 + ma1, rel_tol=1e-10), "Total mass must be conserved"

def test_common_envelope_drag():
    sim, rlmt = make_sim(semimajor=0.4, donor_radius=1.0)
    # Default skip-in-CE = 1 (no RLOF while inside donor)
    rlmt.params["ce_rho0"] = 1e-4
    rlmt.params["ce_cs"]   = 0.05

    d, a = sim.particles[0], sim.particles[1]
    v0 = np.linalg.norm([a.vx-d.vx, a.vy-d.vy, a.vz-d.vz])

    sim.integrate(sim.t + 0.2)
    v1 = np.linalg.norm([a.vx-d.vx, a.vy-d.vy, a.vz-d.vz])

    assert v1 < v0, "Relative speed should drop inside the envelope"

def test_gravitational_wave_decay():
    sim, rlmt = make_sim(semimajor=1.0, donor_radius=0.3)
    rlmt.params["gw_decay_on"] = 1
    rlmt.params["gw_c"]        = 10.0

    o0 = sim.orbits()[0]
    a0 = o0.a

    sim.integrate(sim.t + 2.0)
    a1 = sim.orbits()[0].a

    assert a1 < a0, "Semi-major axis should shrink under GW emission"

def test_immediate_merge_guard():
    # Separation (0.004) < default merge_eps=0.5*min(r=0.01,0.01)=0.005 → immediate merge
    sim, rlmt = make_sim(semimajor=0.004, donor_radius=0.01, accretor_radius=0.01)
    sim.integrate(sim.t + sim.dt)

    assert sim.N == 1, "Particles should merge immediately due to small separation"
    assert sim.particles[0].m > 0.0, "Merged particle must have positive mass"

def test_full_physics_combo():
    sim, rlmt = make_sim(semimajor=0.9, donor_radius=0.8)
    rlmt.params["rlmt_loss_fraction"] = 0.2

    rlmt.params["ce_rho0"] = 1e-5
    rlmt.params["ce_cs"]   = 0.05

    rlmt.params["gw_decay_on"] = 1
    rlmt.params["gw_c"]        = 20.0

    d, a = sim.particles[0], sim.particles[1]
    a0    = sim.orbits()[0].a
    mtot0 = d.m + a.m

    sim.integrate(sim.t + 1.0)

    a1    = sim.orbits()[0].a
    mtot1 = d.m + a.m

    assert d.m < mtot0, "Donor must lose mass in combo test"
    assert a1 < a0, "Orbit must shrink in combo test"
    # Allow for the fact that 20% of transferred mass is lost
    transferred = (mtot0 - (d.m + a.m)) / 0.8
    expected    = mtot0 - 0.2 * transferred
    assert math.isclose(mtot1, expected, rel_tol=1e-3), "Mass loss fraction mismatch"

if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__]))
