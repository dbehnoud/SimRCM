import numpy as np
import cantera as ct
import yaml

class Inputs(object):

    def __init__(self, filename):

        with open(filename, 'r') as f:
            inputs = yaml.safe_load(f)
        
        #vp = inputs['vprofile']
        time_list = [row[0] for row in inputs['vprofile']]
        vol_list = 1e-6*[row[1] for row in inputs['vprofile']]
    
        self.bore = inputs['bore']
        self.t0 = inputs['t0']
        self.z = inputs['znz']
        self.t_wall = inputs['T_wall']
        self.temp0 = inputs['temp0']
        self.p0 = inputs['p0']
        self.v_rcm = 1e-6*vol_list[0]
        self.a_rcm = (np.pi/4)*self.bore**2
        self.mechanism = inputs['mechanism']
        self.mixture = inputs['mixture']
        self.vprofile = {'vproTime': time_list, 'vproVol': vol_list}
        

class VolumeProfile(object):
    """
    Set the velocity of the piston by using a user specified volume
    profile. The initialization and calling of this class are handled
    by the `Func1 <http://cantera.github.io/docs/sphinx/html/cython/zerodim.html#cantera.Func1>`_
    interface of Cantera. Used with the input keyword :ref:`VPRO <VPRO>`
    """

    def __init__(self, keywords, a_rcm):
        """Set the initial values of the arrays from the input keywords.

        The time and volume are read from the input file and stored in
        the ``keywords`` dictionary. The velocity is calculated by
        assuming a unit area and using the forward difference,
        calculated by ``numpy.diff``. This function is only called
        once when the class is initialized at the beginning of a
        problem so it is efficient.

        :param keywords:
            Dictionary of keywords read from the input file
        """

        # The time and volume are stored as lists in the keywords
        # dictionary. The volume is normalized by the first volume
        # element so that a unit area can be used to calculate the
        # velocity.
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
    
    p = np.ones(z)
    p[-1] = p[-1] - bore/2/t0
    all_roots = np.roots(p)
    alpha = all_roots[np.isreal(all_roots)]
    alpha = alpha[alpha>0]
    alpha = np.real(alpha)	
    alpha = alpha*0.99
    
    
    #from numpy.polynomial.polynomial import polyval
    
    class Zone(object):
       
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
            #self.cross_area = np.pi*np.square(self.radius)
                       
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
            qq[x] = k*(r[x].T-t_wall)/(zone[x].thickness)
                
    
    q[1] = -qq[1]
    
    for x in range(2,z+1):	
        q[x] = qq[x-1] - qq[x]
    
    return q

def def_walls(r, env, zone, z, v_factor, keywords, t_wall):
    wq = [0]
    wv = [0]
    q = heat_transfer(z, zone, r, t_wall)  
        
    wq.append(ct.Wall(r[1], env, A = zone[1].surface_area, Q = q[1]))    
    
    for x in range(2,z):	
        wq.append(ct.Wall(r[x], env, A = zone[x].surface_area, Q = q[x]))
        
    wq.append(ct.Wall(r[z], env, A = zone[z].surface_area, Q = q[z]))
    
    for x in range(1,z+1):
        wv.append(ct.Wall(r[x], env, velocity = VolumeProfile(keywords)))
        wv[x].area = v_factor[x] 
    return wq , wv


def modifiy_walls(wq, z, zone):
    
    q = heat_transfer()
    
    for x in range(1,z+1):
        wq[x].set_heat_flux(q[x])
        wq[x].area = zone[x].surface_area
    
    return wq

def find_beta(x, zone):
    
    r = zone[x].radius
    h = zone[x].height  
    
    eqn = np.array([2, (4*r + h), (2*r*h + 2*r**2), (-vsum(x)/np.pi + h*r**2)])
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
        specific_volume_old = r[x].thermo.volume_mass
        gamma = r[x].thermo.cv/r[x].thermo.cp
        
        zone[x].volume = v_old*(r[x].thermo.P/p_rc)**gamma
        r[x].volume = zone[x].volume
        
        specific_volume_new = r[x].volume/v_old * specific_volume_old
        
        v_zones[x] = zone[x].volume
        internal_energy = r[x].thermo.u - (r[x].thermo.P + p_rc)*(r[x].volume-v_old)/r[x].mass/2
        #temperature = r[x].thermo.T*(v_old/zone[x].volume)**(gamma-1)
        contents[x].UV = internal_energy, specific_volume_new
        #contents[x].UV = internal_energy, r[x].thermo.v
        r[x].insert(contents[x])
        
        #accum_v_sum += v_zones[x]
        
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






        
        
    