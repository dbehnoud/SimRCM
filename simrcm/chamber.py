import numpy as np
import cantera as ct
from scipy.optimize import fsolve

f = open('data.txt', 'r')

lines = f.readlines()

val_list = lines[64:-1]

val_dict = {}

time_list = []
vol_list = []

for line in range(len(val_list)):
    str_1 = val_list[line].split('[')[1]
    str_2 = str_1.split(']')[0]
    time, vol = str_2.split(',')

    time_list.append(float(time))
    vol_list.append(float(vol))

    val_dict[float(time)] = float(vol)
	
#for k,v in val_dict.items():
	#print(k,' ',v)
    
class VolumeProfile(object):
    """
    Set the velocity of the piston by using a user specified volume
    profile. The initialization and calling of this class are handled
    by the `Func1 <http://cantera.github.io/docs/sphinx/html/cython/zerodim.html#cantera.Func1>`_
    interface of Cantera. Used with the input keyword :ref:`VPRO <VPRO>`
    """

    def __init__(self, keywords):
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
        self.volume = np.array(keywords['vproVol'])/keywords['vproVol'][0]

        # The velocity is calculated by the forward difference.
        # numpy.diff returns an array one element smaller than the
        # input array, so we append a zero to match the length of the
        # self.time array.
        self.velocity = np.diff(self.volume)/np.diff(self.time)
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
    

bore = 0.0254
t0 = 0.0254/6
alpha = 1
z= 6
v_rcm = val_dict
a_rcm = np.pi/4*bore**2
t=1

keywords = {'vproTime': time_list, 'vproVol': vol_list}


class Zone(object):
    
       
    def __init__(self, i ):
        
        if i < 6:
            from numpy.polynomial.polynomial import polyval
            coef = t0*np.ones(z-i)
              
            self.radius = bore/2-polyval(alpha, coef)
            self.height = v_rcm[0]/a_rcm-2*polyval(alpha, coef)
        else:
            self.radius = bore/2
            self.height = v_rcm[0]/a_rcm
        if i in range(2,z):
            self.thickness = t0*alpha**(z-i)
        else: 
            self.thickness = self.height/2
        self.volume = (np.pi*np.square(self.radius)*self.height)-(np.pi*np.square(self.radius-t0*alpha**(z-i))*(self.height-2*t0*alpha**(z-i)))
        self.outer_surface = 2*np.pi*self.radius*self.height
        
        
gas = ct.Solution('gri30.xml')
gas.TPX = 297.4, 127723, 'H2:0.125,O2:0.0625,N2:0.18125,Ar:0.63125'
gas.transport_model = 'Multi'

env = ct.Reservoir(ct.Solution('air.xml'))


def def_zones(t):
    zone = [0]
    for x in range(1,z+1):
        zone.append(Zone(x))
    return zone

zone = def_zones(t)

contents = [0]
for x in range(1,z+1):
    contents.append(gas)
 
def Reactors(t):
    r=[ct.IdealGasReactor(gas)]
    for x in range(1,z+1):
        r.append(ct.IdealGasReactor(contents[x], volume = zone[x].volume))
    return r

r = Reactors(t)


def walls(t):
    w = [0]
    for x in range(2,z):
        thermal_conductivity = (r[x].thermo.thermal_conductivity + r[x-1].thermo.thermal_conductivity)/2
        u = 2*thermal_conductivity/(zone[x].thickness+zone[x-1].thickness)
        w.append(ct.Wall(r[x], r[x+1], A = zone[x].outer_surface, U = u, velocity = VolumeProfile(keywords)))
    w.append(ct.Wall(r[z], env, velocity = VolumeProfile(keywords)))
    return w

w = walls(t)

def cell_rezoning(t):
    
    accum_v_sum = 0
    pv_zones = [0]
    v_zones = [0]
    
    for x in range(1,z+1):
        pv_zones.append(r[x].thermo.P*r[x].volume)
        v_zones.append(r[x].volume)
    p_rc = sum(pv_zones)/sum(v_zones)
    
    for x in range(1,z+1):
        v_old = r[x].volume
        zone[x].volume = zone[x].volume*(r[x].thermo.P/p_rc)**(r[x].thermo.cv/r[x].thermo.cp)
        r[x].volume = zone[x].volume
        
        internal_energy = r[x].thermo.u - (r[x].thermo.P + p_rc)*(r[x].volume-v_old)/r[x].mass/2
        contents[x].UV = internal_energy, r[x].thermo.v
        
        accum_v_sum += v_zones[x]
        
        def func(beta):
            return accum_v_sum-np.pi*np.square(zone[x].radius + beta)*(zone[x].height + 2*beta)
        
        beta = fsolve(func, 1)
        
        zone[x].radius = zone[x].radius + beta
        zone[x].height = zone[x].height + 2*beta
        
netw = ct.ReactorNet(r)

time = []
temperature = []
pressure = []
volume = []
mass_fractions = []

while netw.time < 0.05:
    time.append(netw.time)
    for j in range(1,z+1):
        temperature.append(r[j].thermo.T)
        pressure.append(r[j].thermo.P)
        volume.append(r[j].volume)
        mass_fractions.append(r[j].Y)
    netw.step()
    cell_rezoning(t)

        
        
        
        

        
    
           


    
    
        
        
        
        
        
        





