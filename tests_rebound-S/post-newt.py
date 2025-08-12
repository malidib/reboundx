import numpy as np
import rebound
import reboundx

# ---------- Units ----------
# Astronomical units: length=AU, time=yr, mass=Msun
G = 4*np.pi**2                    # AU^3 / (Msun * yr^2)
c_au_per_yr = 63241.077           # speed of light in AU/yr

# ---------- Simulation ----------
sim = rebound.Simulation()
sim.units = ('AU','yr','Msun')    # purely for your own bookkeeping/logging
sim.integrator = "ias15"          # velocity-dependent forces -> IAS15 recommended
sim.G = G

# Equal-mass compact binary with mild eccentricity
m1 = 30.0       # Msun
m2 = 30.0
a0 = 0.01       # AU
e0 = 0.10
sim.add(m=m1, x=-a0*m2/(m1+m2), y=0, z=0, vx=0, 
        vy= np.sqrt(G*(m1+m2)/a0)*np.sqrt((1-e0)/(1+e0))*m2/(m1+m2),
        vz=0)
sim.add(m=m2, x= a0*m1/(m1+m2), y=0, z=0, vx=0, 
        vy=-np.sqrt(G*(m1+m2)/a0)*np.sqrt((1-e0)/(1+e0))*m1/(m1+m2),
        vz=0)

# Move to center of mass (always good practice)
sim.move_to_com()

# ---------- REBOUNDx: attach your PN force ----------
rebx = reboundx.Extras(sim)

# IMPORTANT: your C force must be compiled/installed so that REBOUNDx can find it by name.
# For example, place post_newtonian.c in your local reboundx/effects build and reinstall.
# Now load it like any other force:
pn = rebx.load_force("post_newtonian")     # <-- your module name
rebx.add_force(pn)

# Required parameter: speed of light in *your* units
pn.params["c"] = c_au_per_yr              # AU/yr

# Optional switches (default 1). You can toggle pieces as needed.
pn.params["pn_15PN"] = 1
pn.params["pn_2PN"]  = 1
pn.params["pn_25PN"] = 1

# ---------- (Optional) Spins ----------
# Your C force expects a particle param "pn_spin" with dimensions mass*L^2/T.
# Here we give modest, misaligned spins to illustrate precession.
S1 = np.array([0.0, 0.0, 0.05])*m1*a0**2/1.0    # toy numbers, pick your model
S2 = np.array([0.0, 0.0, -0.05])*m2*a0**2/1.0
sim.particles[0].params["pn_spin"] = S1
sim.particles[1].params["pn_spin"] = S2

# ---------- Integrate and record orbital decay ----------
def rel_orbit(sim):
    p0, p1 = sim.particles[0], sim.particles[1]
    orb = p1.orbit(primary=p0)        # current API
    return orb.a, orb.e, orb.omega

times = np.linspace(0, 3.0, 601)  # 3 years should show visible inspiral at a0=0.01 AU
a_hist, e_hist, w_hist = [], [], []

for t in times:
    sim.integrate(t)
    a, e, w = rel_orbit(sim)
    a_hist.append(a); e_hist.append(e); w_hist.append(w)

print(f"Initial a={a_hist[0]:.6f} AU -> final a={a_hist[-1]:.6f} AU")
print(f"Initial e={e_hist[0]:.4f}  -> final e={e_hist[-1]:.4f}")
print(f"Pericenter advance Δω ≈ {(w_hist[-1]-w_hist[0]):.4f} rad")

# (Optional) quick sanity check for monotonic inspiral with 2.5PN on:
assert a_hist[-1] < a_hist[0], "a did not decrease; check pn_25PN and c units."
