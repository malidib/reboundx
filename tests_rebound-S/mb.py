import numpy as np
import rebound, reboundx

# ---------- Units ----------
sim = rebound.Simulation()
sim.units = ('AU','yr','Msun')
sim.integrator = "ias15"
sim.G = 4*np.pi**2
rebx = reboundx.Extras(sim)

# ---------- Binary ICs ----------
AU_per_Rsun = 1.0/215.032
m1 = m2 = 1.0
a0, e0 = 0.05, 0.05
R1 = R2 = 1.0*AU_per_Rsun

sim.add(m=m1, r=R1, x=-a0*m2/(m1+m2),
        vy=+np.sqrt(sim.G*(m1+m2)/a0)*np.sqrt((1-e0)/(1+e0))*m2/(m1+m2))
sim.add(m=m2, r=R2, x=+a0*m1/(m1+m2),
        vy=-np.sqrt(sim.G*(m1+m2)/a0)*np.sqrt((1-e0)/(1+e0))*m1/(m1+m2))
sim.move_to_com()

# ---------- Spins & structure ----------
k2_moi = 0.07
I1 = k2_moi*m1*R1**2
I2 = k2_moi*m2*R2**2
P_day = 1.0/365.25
omega0 = 2*np.pi/P_day

for i,(I,R) in enumerate([(I1,R1),(I2,R2)]):
    p = sim.particles[i]
    p.params["Omega"] = np.array([0.0,0.0,omega0])
    p.params["I"]     = I
    p.params["k2"]    = 0.02           # Love number
    p.params["tau"]   = 0.05/365.25    # only used if tides_spin is available

    # magnetic_braking inputs
    p.params["mb_on"]         = 1
    p.params["mb_convective"] = 1
    p.params["mb_tau_conv"]   = 12.5/365.25

# ---------- Add magnetic_braking ----------
mb = rebx.load_operator("magnetic_braking")
rebx.add_operator(mb)
mb.params["mb_Msun"] = 1.0
mb.params["mb_Rsun"] = AU_per_Rsun
mb.params["mb_year"] = 1.0
mb.params["mb_Rossby_sat"] = 0.1
# mb.params["mb_K"] = 1e49  # (optional) stronger for fast demos

# ---------- Try tides_spin; fallback to tides_constant_time_lag ----------
try:
    tides = rebx.load_force('tides_constant_time_lag')
    rebx.add_force(tides)
    use_spin = True
    # spins will evolve self-consistently; nothing else to do
except Exception:
    tides = rebx.load_operator("tides_constant_time_lag")
    rebx.add_operator(tides)
    use_spin = False
    # set tctl params (older operator; no spin evolution)
    for i,p in enumerate(sim.particles[:2]):
        p.params["tctl_k2"]  = 0.02
        p.params["tctl_tau"] = 0.05/365.25
        # OmegaMag is read by tctl; we will sync it from our Omega vector each step
        Om = p.params["Omega"]; p.params["OmegaMag"] = np.linalg.norm(Om)

def elems():
    p1=sim.particles[1]
    p0=sim.particles[0]
    o = p1.orbit(primary=p0)
    return o.a, o.e

def sync_spin_to_tctl():
    if use_spin: return
    for p in sim.particles[:2]:
        Om = p.params["Omega"]
        p.params["OmegaMag"] = float(np.linalg.norm(Om))

def prot_days(i):
    w = np.linalg.norm(sim.particles[i].params["Omega"])
    return (2*np.pi/w)*365.25

# ---------- Integrate ----------
t_end = 2e3
N = 10
times = np.linspace(0.0, t_end, N)
a_hist, e_hist = [], []
P1_hist, P2_hist = [], []

for t in times:
    sim.integrate(t)
    sync_spin_to_tctl()  # keep tctl aware of magnetic braking spin-down
    a,e = elems()
    a_hist.append(a); e_hist.append(e)
    P1_hist.append(prot_days(0)); P2_hist.append(prot_days(1))

print(("Using tides_spin" if use_spin else "Using tides_constant_time_lag (with Omega→OmegaMag sync)"))
print(f"a: {a_hist[0]:.6f} AU  ->  {a_hist[-1]:.6f} AU")
print(f"e: {e_hist[0]:.4f}     ->  {e_hist[-1]:.4f}")
print(f"P1: {P1_hist[0]:.2f} d ->  {P1_hist[-1]:.2f} d | P2: {P2_hist[0]:.2f} d -> {P2_hist[-1]:.2f} d")
