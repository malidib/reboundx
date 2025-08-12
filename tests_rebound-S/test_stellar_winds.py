import rebound
import reboundx

sim = rebound.Simulation()
sim.units = ('AU', 'yr', 'Msun')
sim.integrator='whfast'
sim.dt=5000
sim.add(m=100.0)
star = sim.particles[0]


rebx = reboundx.Extras(sim)
#swml = rebx.load_operator('stellar_wind_mass_loss')
#rebx.add_operator(swml)

sse = rebx.load_operator('stellar_evolution_sse')
rebx.add_operator(sse)
sim.integrate(sim.t + 1)
print (star.m,star.r, star.params["sse_L"])
tdw = rebx.load_operator('thermally_driven_winds')
rebx.add_operator(tdw)
print (star.m,star.r, star.params["sse_L"])
star.params['tdw_eta'] = 100000
tdw.params['tdw_const'] = 1e1000

sim.integrate(sim.t + 2.0e8)
#star.params['swml_L']=1
print (star.m,star.r, star.params["sse_L"])

#star.params['swml_eta'] = 1



#star.params['tdw_T'] = 4e7
#star.params['tdw_R']= star.r
#star.params["sse_R_coeff"] = 10.5
#star.params['sse_L_coeff']=1.
#rebx.set_operator_parameter(op, "sse_R_exp", 0.85)
#rebx.set_operator_parameter(op, "sse_L_coeff", 1.2)
#print (star.m,star.r, star.params["swml_L"])


#star.params['tdw_T'] = 1.5e6
#star.params['tdw_R'] = star.r

initial_mass = star.m
for _ in range(5):
        sim.integrate(sim.t + 2.0e8)



print (star.m,star.r, star.params["sse_L"])
#print (star.params["swml_R"])