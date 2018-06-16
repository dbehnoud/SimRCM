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
    
    assert zone[1].radius == 9.27566302
    assert zone[1].height == 14.55132604
    assert zone[1].thickness == 9.27566302
    assert zone[1].volume == 3933.16684353
    assert zone[1].surface_area == 1388.65382985
    
    assert zone[3].radius == 17.2831974
    assert zone[3].height == 30.56639479
    assert zone[3].thickness == 2.94741118
    assert zone[3].volume == 12755.12408476
    assert zone[3].surface_area == 5196.15621734
    
    assert zone[5].radius == 20.0
    assert zone[5].height == 36.0
    assert zone[5].thickness == 1.0
    assert zone[5].volume == 6679.02598153
    assert zone[5].surface_area == 7037.167544041136
    
