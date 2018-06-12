import sim_tools 
import cantera as ct

inputs = sim_tools.Inputs('assets\example_input.yaml')

bore = inputs.bore
t0 = inputs.t0
z = inputs.z
v_rcm = inputs.v_rcm
a_rcm = inputs.a_rcm
t_wall = inputs.T_wall

zone = sim_tools.def_zones()

v_factor = [0]
for x in range(1,z+1):
    v_factor.append(zone[x].r_volume*a_rcm)
    
r = sim_tools.def_reactors()
wq , wv = sim_tools.def_walls()

netw =[0]        
for x in range(1,z+1):        
    netw.append(ct.ReactorNet([r[x]]))
    
#while netw[x].time < 0.05:
    for x in range(1,z+1):
        netw[x].step()
    zone, r = sim_tools.cell_rezone() 
    wq = sim_tools.modify_walls(wq) 