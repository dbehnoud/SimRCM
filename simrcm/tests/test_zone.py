import numpy as np
import simrcm.sim_tools
from numpy.polynomial.polynomial import polyval
 
    
def test_zone():
    t0 = 1
    z =5
    bore = 40
    alpha = 2
    a_rcm = np.pi/4*bore**2
    v_rcm = a_rcm*36
    
    zone = simrcm.sim_tools.def_zones(z, bore, t0, v_rcm, a_rcm)
    
    assert zone[1].radius == 5
    assert zone[1].height == 6
    assert zone[1].thickness == 5
    assert zone[1].volume == np.pi*25*6
    assert zone[1].surface_area == 2*np.pi*25 + 2*np.pi*30
    
    assert zone[3].radius == 17
    assert zone[3].height == 30
    assert zone[3].thickness == 4
    assert zone[3].volume == np.pi*17**2*30 - np.pi*13**2*22
    assert zone[3].surface_area == 2*np.pi*17**2 + 2*np.pi*17*30
    
    assert zone[5].radius == 20
    assert zone[5].height == 36
    assert zone[5].thickness == 1
    assert zone[5].volume == np.pi*400*36 - np.pi*19**2*34
    assert zone[5].surface_area == 2*np.pi*400 + 2*np.pi*20*36
    
