import numpy as np

def def_zones(t0,z,bore,alpha,a_rcm,v_total):
    
    from numpy.polynomial.polynomial import polyval
    
    class Zone(object):
       
        def __init__(self, i ):
            
            if i < z:
                coef = t0*np.ones(z-i)
                  
                self.radius = bore/2-polyval(alpha, coef)
                self.height = v_total/a_rcm-2*polyval(alpha, coef)
            else:
                self.radius = bore/2
                self.height = v_total/a_rcm     
        
            if i == 1:
                self.thickness = self.radius
            else: 
                self.thickness = t0*alpha**(z-i)
                
            if i == 1:
                self.volume = np.pi*np.square(self.radius)*self.height
            else:
                self.volume = (np.pi*np.square(self.radius)*self.height)-(np.pi*np.square(self.radius-self.thickness)*(self.height-2*self.thickness))
            
            self.surface_area = 2*np.pi*self.radius*self.height + 2*np.pi*np.square(self.radius)
            self.cross_area = np.pi*np.square(self.radius)
         
    zone = [0]
    for x in range(1,z+1):
        zone.append(Zone(x))
    return zone


def zone_test():
    
    t0 = 1
    z =5
    bore = 40
    alpha = 2
    a_rcm = np.pi/4*bore**2
    v_total = a_rcm*36
    
    zone = def_zones(t0,z,bore,alpha,a_rcm,v_total)
    
    assert zone[1].radius == 5
    assert zone[1].height == 6
    assert zone[1].thickness == 5
    assert zone[1].volume == np.pi*25*6
    assert zone[1].surface_area == 2*np.pi*25 + 2*np.pi*30
    assert zone[1].cross_area == np.pi*25
    
    assert zone[3].radius == 17
    assert zone[3].height == 30
    assert zone[3].thickness == 4
    assert zone[3].volume == np.pi*17**2*30 - np.pi*13**2*22
    assert zone[3].surface_area == 2*np.pi*17**2 + 2*np.pi*17*30
    assert zone[3].cross_area == np.pi*17**2
    
    assert zone[5].radius == 20
    assert zone[5].height == 36
    assert zone[5].thickness == 1
    assert zone[5].volume == np.pi*400*36 - np.pi*19**2*34
    assert zone[5].surface_area == 2*np.pi*400 + 2*np.pi*20*36
    assert zone[5].cross_area == np.pi*400