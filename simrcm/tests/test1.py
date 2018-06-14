import yaml
import numpy

class Inputs(object):

    def __init__(self, filename):

        with open(filename, 'r') as f:
            inputs = yaml.safe_load(f)
        
        #vp = inputs['vprofile']
        time_list = [row[0] for row in inputs['vprofile']]
        vol_list = [row[1] for row in inputs['vprofile']]
    
        self.bore = inputs['bore']
        self.t0 = inputs['t0']
        self.z = inputs['znz']
        self.t_wall = inputs['T_wall']
        self.v_rcm = vol_list[0]
        self.a_rcm = (numpy.pi/4)*self.bore**2
        self.mechanism = inputs['mechanism']
        self.vprofile = {'vproTime': time_list, 'vproVol': vol_list}
        
inputs = Inputs('example_input.yaml')

assert inputs.z == 10
assert inputs.bore == 0.0508 
assert inputs.v_rcm == 5.47669375000E+002         
assert inputs.vprofile['vproVol'][2] == 5.43427034574E+002     
        