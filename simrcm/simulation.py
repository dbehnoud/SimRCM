import simrcm.sim_tools 
import cantera as ct
import numpy as np

def simulation1(filename):
    """
    The main algorithim for running a simulation case and evaluationg the ignition delay time
    
    Returns
    -------
    ignition_delay : Igition delay time
    pressure: List of reactor pressure at each time step
    temperature: List of reactor temperature at each time step
    time: List of time steps used in simulation
    """

    inputs = simrcm.sim_tools.Inputs(filename)
    
    bore = inputs.bore
    t0 = inputs.t0
    z = inputs.z
    v_rcm = inputs.v_rcm
    a_rcm = inputs.a_rcm
    t_wall = inputs.t_wall
    temp0 = inputs.temp0
    p0 = inputs.p0
    mechanism = inputs.mechanism
    mixture = inputs.mixture
    keywords = inputs.vprofile
    
    # Zone definition
    zone = simrcm.sim_tools.def_zones(z, bore, t0, v_rcm, a_rcm)
    # Storing the velocity factor 
    v_factor = [0]
    for x in range(1,z+1):
        v_factor.append(zone[x].r_volume*a_rcm)
    # Reactors definition  
    r, env, contents = simrcm.sim_tools.def_reactors(z, zone, temp0, p0, mechanism, mixture)
    # Walss definition
    wq , wv = simrcm.sim_tools.def_walls(r, env, zone, z, v_factor, keywords, t_wall, a_rcm)
    # Integration
    netw =[0]        
    for x in range(1,z+1):        
        netw.append(ct.ReactorNet([r[x]]))
        
    time = []
    pressure = []
    temperature = []
            
    while netw[1].time < 0.05:
        time.append(netw[1].time)
        temperature.append(r[1].thermo.T)
        pressure.append(r[1].thermo.P)
        for x in range(1,z+1):
            netw[x].step()   
        zone, r = simrcm.sim_tools.cell_rezone(z, r, zone, contents) 
        wq = simrcm.sim_tools.modify_walls(wq, z, zone, r, t_wall) 
    
    # Ignition delay time evaluation
    dpdt = np.gradient(pressure, time, edge_order=2) 
    ignition_delay = time[np.argmax(dpdt)]
    
    
    
    return ignition_delay, pressure, temperature, time
