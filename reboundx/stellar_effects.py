"""Simplified stellar physics utilities for REBOUNDx.
These Python-level implementations provide approximate versions of
common stellar effects and are intended for demonstration only."""

import math
import rebound
import reboundx

__all__ = [
    'add_roche_lobe_mass_transfer',
    'add_tidal_dissipation',
    'add_gravitational_wave_decay',
    'add_supernova_event',
    'add_magnetic_braking',
    'add_common_envelope_drag',
    'add_stellar_wind_mass_loss',
    'add_simple_stellar_evolution',
    'add_non_conservative_mass_loss',
    'add_gravitational_harmonics_evolving',
]


def add_roche_lobe_mass_transfer(rebx, donor, accretor, mdot, beta=1.0):
    """Apply a simple Roche-lobe mass transfer between two bodies."""
    op = rebx.create_operator('roche_lobe_transfer')
    op.operator_type = 'updater'
    donor_i = donor
    acc_i = accretor
    mdot_val = mdot
    beta_val = beta

    def step(sim, operator, dt, donor_i=donor_i, acc_i=acc_i, mdot_val=mdot_val, beta_val=beta_val):
        d = donor_i
        a = acc_i
        md = mdot_val * dt
        b = beta_val
        p_d = sim.contents.particles[d]
        p_a = sim.contents.particles[a]
        p_d.m -= md
        p_a.m += b * md
    op.step_function = step
    rebx.add_operator(op)
    return op


def add_tidal_dissipation(rebx):
    """Load the built-in constant time-lag tidal forces with spin coupling."""
    f1 = rebx.load_force('tides_constant_time_lag')
    f2 = rebx.load_force('tides_spin')
    rebx.add_force(f1)
    rebx.add_force(f2)
    return f1, f2


def add_gravitational_wave_decay(rebx, body1, body2, c=299792458.0):
    """Add a crude gravitational-wave orbital decay operator."""
    op = rebx.create_operator('gw_decay')
    op.operator_type = 'updater'
    b1 = body1
    b2 = body2
    cval = c

    def step(sim, operator, dt, b1=b1, b2=b2, cval=cval):
        i = b1
        j = b2
        c = cval
        p1 = sim.contents.particles[i]
        p2 = sim.contents.particles[j]
        orb = rebound.Orbit.from_binary(sim, p1, p2)
        mu = p1.m + p2.m
        da_dt = -(64./5.) * p1.m * p2.m * mu / (c**5 * orb.a**3 * (1-orb.e**2)**3.5)
        orb.a += da_dt * dt
        newp1, newp2 = orb.to_binary(sim.G, p1.m, p2.m)
        sim.contents.particles[i] = newp1
        sim.contents.particles[j] = newp2
    op.step_function = step
    rebx.add_operator(op)
    return op


def add_supernova_event(rebx, body, t_event, m_final, vkick=(0,0,0)):
    """Trigger instantaneous mass loss and kick at t_event."""
    op = rebx.create_operator('supernova_event')
    op.operator_type = 'updater'
    b = body
    tev = t_event
    mf = m_final
    kx, ky, kz = vkick
    done = False

    def step(sim, operator, dt):
        nonlocal done
        if done:
            return
        if sim.contents.t >= tev:
            idx = int(b)
            p = sim.contents.particles[idx]
            p.m = mf
            p.vx += kx
            p.vy += ky
            p.vz += kz
            done = True
    op.step_function = step
    rebx.add_operator(op)
    return op


def add_magnetic_braking(rebx, body, K):
    """Apply a simple magnetic braking torque on a star."""
    op = rebx.create_operator('magnetic_braking')
    op.operator_type = 'updater'
    b = body
    Kval = K

    def step(sim, operator, dt, b=b, Kval=Kval):
        idx = b
        p = sim.contents.particles[idx]
        Omega = p.params.get('Omega', 0.0)
        Omega += -Kval * Omega**3 * dt
        p.params['Omega'] = Omega
    op.step_function = step
    rebx.add_operator(op)
    return op


def add_common_envelope_drag(rebx, body, Cd):
    """Simple dynamical friction inside a common envelope."""
    force = rebx.create_force('common_envelope_drag')
    force.force_type = 'vel'
    force.body = body
    force.Cd = Cd

    def f(sim, force, particles, N):
        idx = int(force.body)
        Cd = force.Cd
        p = particles[idx]
        p.ax += -Cd * p.vx
        p.ay += -Cd * p.vy
        p.az += -Cd * p.vz
    force.update_accelerations = f
    rebx.add_force(force)
    return force


def add_stellar_wind_mass_loss(rebx, body, mdot):
    """Continuous isotropic stellar wind mass loss."""
    op = rebx.create_operator('stellar_wind_mass_loss')
    op.operator_type = 'updater'
    b = body
    md = mdot

    def step(sim, operator, dt, b=b, md=md):
        idx = b
        p = sim.contents.particles[idx]
        p.m -= md * dt
    op.step_function = step
    rebx.add_operator(op)
    return op


def add_simple_stellar_evolution(rebx, body, mass0, mass_final, t_end):
    """Linear decrease in stellar mass until t_end."""
    op = rebx.create_operator('simple_stellar_evolution')
    op.operator_type = 'updater'
    b = body
    m0 = mass0
    mf = mass_final
    tend = t_end

    def step(sim, operator, dt):
        idx = b
        p = sim.contents.particles[idx]
        t = sim.contents.t
        if t < tend:
            p.m = m0 + (mf - m0) * t / tend
        else:
            p.m = mf
    op.step_function = step
    rebx.add_operator(op)
    return op


def add_non_conservative_mass_loss(rebx, donor, accretor, mdot, beta):
    """Mass transfer where only a fraction beta is accreted."""
    op = rebx.create_operator('non_conservative_mass_loss')
    op.operator_type = 'updater'
    donor_i = donor
    acc_i = accretor
    mdot_val = mdot
    beta_val = beta

    def step(sim, operator, dt):
        d = donor_i
        a = acc_i
        md = mdot_val * dt
        b = beta_val
        p_d = sim.contents.particles[d]
        p_a = sim.contents.particles[a]
        p_d.m -= md
        p_a.m += b * md
    op.step_function = step
    rebx.add_operator(op)
    return op


def add_gravitational_harmonics_evolving(rebx, body, J2, J4=0.0):
    """Wrapper for the built-in gravitational_harmonics force."""
    force = rebx.load_force('gravitational_harmonics')
    force.params['J2'] = J2
    force.params['J4'] = J4
    force.params['central'] = body
    rebx.add_force(force)
    return force

