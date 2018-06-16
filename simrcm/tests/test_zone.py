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
    
    assert zone[1].radius == pytest.approx(9.27566302, rel=1e-1)
    assert zone[1].height == pytest.approx(14.55132604, rel=1e-1)
    assert zone[1].thickness == pytest.approx(9.27566302, rel=1e-1)
    assert zone[1].volume == pytest.approx(3933.16684353, rel=1e-1)
    assert zone[1].surface_area == pytest.approx(1388.65382985, rel=1e-1)
    
    assert zone[3].radius == pytest.approx(17.2831974, rel=1e-1)
    assert zone[3].height == pytest.approx(30.56639479, rel=1e-1)
    assert zone[3].thickness == pytest.approx(2.94741118, rel=1e-1)
    assert zone[3].volume == pytest.approx(12755.12408476, rel=1e-1)
    assert zone[3].surface_area == pytest.approx(5196.15621734, rel=1e-1)
    
    assert zone[5].radius == pytest.approx(20.0, rel=1e-1)
    assert zone[5].height == pytest.approx(36.0, rel=1e-1)
    assert zone[5].thickness == pytest.approx(1.0, rel=1e-1)
    assert zone[5].volume == pytest.approx(6679.02598153, rel=1e-1)
    assert zone[5].surface_area == pytest.approx(7037.167544041136, rel=1e-1)
    
