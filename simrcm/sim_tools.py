import numpy as np
import cantera as ct
import yaml

class Inputs(object):
    """
    Reads the YAML-encoded input file and sets the input values
    
    Parameters
    ----------
    filename : String filename of the input file
        
    Attributes
    ----------
    bore : Bore of the reaction chamber
    t0 : Thickness of outermost zone
    z : Number of zones
    t_wall : Wall temperature
    temp0 : The initial temperature of the simulation
    p0 : The initial pressure of the simulation
    v_rcm : Volume of the reactor at which the compression starts
    a_rcm : Cross sectional area of the reaction chamber
    mechanism : Chemical kinetic mechanism
    mixture : Initial chemical composition
    vprofile : The volume trace to be used for the simulation
        
    """

    def __init__(self, filename):
        
        with open(filename, 'r') as f:
            inputs = yaml.safe_load(f)
        
        time_list = [row[0] for row in inputs['vprofile']]
        vol_list = []
        for x in [row[1] for row in inputs['vprofile']]:
            vol_list.append(1e-6*float(x))    
              
        self.bore = inputs['bore']
        self.t0 = inputs['t0']
        self.z = inputs['znz']
        self.t_wall = inputs['T_wall']
        self.temp0 = inputs['temp0']
        self.p0 = inputs['p0']
        self.v_rcm = vol_list[0]
        self.a_rcm = (np.pi/4)*self.bore**2
        self.mechanism = inputs['mechanism']
        self.mixture = inputs['mixture']
        self.vprofile = {'vproTime': time_list, 'vproVol': vol_list}
        

class VolumeProfile(object):
    """
    Set the velocity of the piston by using a user specified volume
    profile.
    The time and volume are read from the input file and stored in
        the ``keywords`` dictionary. 

        :param keywords:
            Dictionary of keywords read from the input file
    """

    def __init__(self, keywords, a_rcm):
        
        # The time and volume are stored as lists in the keywords
        # dictionary. 
        self.time = np.array(keywords['vproTime'])
        self.volume = np.array(keywords['vproVol'])

        # The velocity is calculated by the forward difference.
        # numpy.diff returns an array one element smaller than the
        # input array, so we append a zero to match the length of the
        # self.time array.
        self.velocity = np.diff(self.volume/a_rcm)/np.diff(self.time)
        self.velocity = np.append(self.velocity, 0)

    def __call__(self, t):
        """Return the velocity when called during a time step.

        Using linear interpolation, determine the velocity at a given
        input time ``t``.

        :param t:
        Input float, current simulation time.
        """

        if t < self.time[-1]:
             # prev_time_point is the previous value in the time array
             # after the current simulation time
             prev_time_point = self.time[self.time <= t][-1]
             # index is the index of the time array where
             # prev_time_point occurs
             index = np.where(self.time == prev_time_point)[0][0]
             
             return self.velocity[index]
        else:
             return 0
         

def def_zones(z, bore, t0, v_rcm, a_rcm):
    """
    Generates the initial mesh and defines zones
    
    Parameters
    ----------
    z : Number of zones
    bore : Bore of the reaction chamber
    t0 : Thickness of outermost zone
    v_rcm : Volume of the reactor at which the compression starts
    a_rcm : Cross sectional area of the reaction chamber
    
    Returns
    -------
    zone : List of Zone objects
    """
    # Root finding algorithm for finding the growth factor (alpha)
    p = np.ones(z)
    p[-1] = p[-1] - bore/2/t0
    all_roots = np.roots(p)
    alpha = all_roots[np.isreal(all_roots)]
    alpha = alpha[alpha>0]
    alpha = np.real(alpha)	
    alpha = alpha*0.99
     
    
    class Zone(object):
        """
        stores the geometrical information of zones
        
        Parameters
        ----------
        i : index of a single zone
            
        Attributes
        ----------
        radius : Radius of a single zone
        height : Thickness of a single zone
        volume : Volume of a single zone
        thickness : Thickness of a single zone
        surface_area : Outer surface area of a single zone
        """
       
        def __init__(self, i ):
            
            if i < z:
                coef = t0*np.ones(z-i)
                  
                self.radius = bore/2-np.polyval(coef, alpha)
                self.height = v_rcm/a_rcm-2*np.polyval(coef, alpha)
            else:
                self.radius = bore/2
                self.height = v_rcm/a_rcm     
        
            if i == 1:
                self.thickness = self.radius
            else: 
                self.thickness = t0*alpha**(z-i)
                
            if i == 1:
                self.volume = np.pi*np.square(self.radius)*self.height
            else:
                self.volume = (np.pi*np.square(self.radius)*self.height)-(np.pi*np.square(self.radius-self.thickness)*(self.height-2*self.thickness))
            
            self.r_volume = self.volume/v_rcm
            
            self.surface_area = 2*np.pi*self.radius*self.height + 2*np.pi*np.square(self.radius)
            
                       
    zone = [0]
    for x in range(1,z+1):
        zone.append(Zone(x))
        
    return zone


def def_reactors(z, zone, temp0, p0, mechanism, mixture):
    
    gas = ct.Solution(mechanism)
    gas.TPX = temp0, p0, mixture
    gas.transport_model = 'Multi'

    env = ct.Reservoir(ct.Solution(mechanism))
    
    contents = [0]
    for x in range(1,z+1):
        contents.append(gas)
    
    r=[ct.IdealGasReactor(gas)]
    for x in range(1,z+1):
        r.append(ct.IdealGasReactor(contents[x], volume = zone[x].volume))
          
    return r, env, contents

def heat_transfer(z, zone, r, t_wall):
    qq = np.zeros(z+1)
    q = np.zeros(z+1)
    
    for x in range(1,z+1):
        
        if x < z:
            k = (r[x].thermo.thermal_conductivity + r[x+1].thermo.thermal_conductivity)/2
            qq[x] = k*(r[x].T-r[x+1].T)/((zone[x].thickness+zone[x+1].thickness)/2) 
        else:
            k = (r[x].thermo.thermal_conductivity + r[0].thermo.thermal_conductivity)/2
            qq[x] = 2*k*(r[x].T-t_wall)/(zone[x].thickness)
                
    
    q[1] = -qq[1]
    
    for x in range(2,z+1):	
        q[x] = qq[x-1] - qq[x]
    
    return q

def def_walls(r, env, zone, z, v_factor, keywords, t_wall, a_rcm):
    wq = [0]
    wv = [0]
    q = heat_transfer(z, zone, r, t_wall)  
        
    wq.append(ct.Wall(r[1], env, A = zone[1].surface_area, Q = q[1]))    
    
    for x in range(2,z):	
        wq.append(ct.Wall(r[x], env, A = zone[x].surface_area, Q = q[x]))
        
    wq.append(ct.Wall(r[z], env, A = zone[z].surface_area, Q = q[z]))
    
    for x in range(1,z+1):
        wv.append(ct.Wall(r[x], env, velocity = VolumeProfile(keywords,a_rcm)))
        wv[x].area = v_factor[x] 
    return wq , wv


def modify_walls(wq, z, zone, r, t_wall):
    
    q = heat_transfer(z, zone, r, t_wall)
    
    for x in range(1,z+1):
        wq[x].set_heat_flux(q[x])
        wq[x].area = zone[x].surface_area
    
    return wq

def find_beta(x, zone):
    
    r = zone[x].radius
    h = zone[x].height  
    
    eqn = np.array([2, (4*r + h), (2*r*h + 2*r**2), (-vsum(x, zone)/np.pi + h*r**2)])
    eqn_roots = np.roots(eqn)
    beta_it = eqn_roots[np.isreal(eqn_roots)]
    beta_it = np.real(beta_it)
    beta_it = beta_it[beta_it>-r]
    beta_it = beta_it[beta_it<r]
    if len(beta_it) == 1:
        beta = float(beta_it)
    else:
        beta = 0
        
    return beta

def cell_rezone(z, r, zone, contents):
    
    #accum_v_sum = 0
    pv_zones = [0]
    v_zones = [0]
    
    for x in range(1,z+1):
        pv_zones.append(r[x].thermo.P*r[x].volume)
        v_zones.append(r[x].volume)
    p_rc = sum(pv_zones)/sum(v_zones)
    
    for x in range(1,z+1):
        v_old = r[x].volume
        x_old = r[x].thermo.X
        specific_volume_old = r[x].thermo.volume_mass
        gamma = r[x].thermo.cv_mole/r[x].thermo.cp_mole
        
        zone[x].volume = v_old*(r[x].thermo.P/p_rc)**gamma
        r[x].volume = zone[x].volume
        
        specific_volume_new = r[x].volume/v_old * specific_volume_old
        
        v_zones[x] = zone[x].volume
        #internal_energy = r[x].thermo.u - (r[x].thermo.P + p_rc)*(r[x].volume-v_old)/r[x].mass/2
        temperature = r[x].thermo.T*(v_old/zone[x].volume)**(gamma-1)
        
        #contents[x].UV = internal_energy, specific_volume_new
        contents[x].TPX = temperature, p_rc, x_old
        #contents[x].UV = r[x].thermo.u, specific_volume_new
        #r[x].insert(contents[x])
        
               
        beta = find_beta(x, zone)
        
        #print("beta=%s" % beta)
        
        zone[x].radius += beta
        zone[x].height += 2*beta
        
        if x == 1:
            zone[x].thickness = zone[x].radius	
        else:
            zone[x].thickness = zone[x].radius-zone[x-1].radius
            
        return zone, r
    
def vsum(x, zone):
    v = 0.0
    for x in range(1,x+1):
        v += zone[x].volume
    return v


class State(object):

    def __init__(self):
        
        self.pressure = []
        self.temperature = []
        
def get_species_mass(i, r, z):
        oh = 0; rm = 0
        for x in range(1,z+1):
            oh += r[x].thermo.Y[i]*r[x].mass
            rm += r[x].mass
        return oh/rm
                    
                    
                    
                    
                    
                    
    






        
        
    