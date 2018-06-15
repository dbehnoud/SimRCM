import sim_tools 
import cantera as ct
import numpy as np

def simulation1(file):

    inputs = sim_tools.Inputs(file)
    
    bore = inputs.bore
    t0 = inputs.t0
    z = inputs.z
    v_rcm = inputs.v_rcm
    a_rcm = inputs.a_rcm
    t_wall = inputs.T_wall
    temp0 = inputs.temp0
    p0 = inputs.p0
    mechanism = inputs.mechanism
    mixture = inputs.mixture
    keywords = inputs.vprofile
    
    zone = sim_tools.def_zones(z, bore, t0, v_rcm, a_rcm)
    
    v_factor = [0]
    for x in range(1,z+1):
        v_factor.append(zone[x].r_volume*a_rcm)
        
    r, env, contents = sim_tools.def_reactors(z, zone, temp0, p0, mechanism, mixture)
    wq , wv = sim_tools.def_walls(r, env, zone, z, v_factor, keywords, t_wall)
    
    netw =[0]        
    for x in range(1,z+1):        
        netw.append(ct.ReactorNet([r[x]]))
        
    time = []
    pressure = []
    temperature = []
        
    while netw[x].time < 0.05:
        temperature.append(r[1].thermo.T)
        pressure.append(r[1].thermo.P)
        for x in range(1,z+1):
            netw[x].step()   
        zone, r = sim_tools.cell_rezone(z, r, zone, contents) 
        wq = sim_tools.modify_walls(wq, z, zone) 
        
    dpdt = np.gradient(pressure, time, edge_order=2) 
    ignition_delay = time[np.argmax(dpdt)]
    
    import matplotlib.pyplot as plt
    
    plt.figure()
    plt.plot(time, pressure)
    plt.ylabel('Pressure [Pa]')
    plt.xlabel('Time [s]')
    plt.grid(True, which='both', axis='x');
    
    return ignition_delay