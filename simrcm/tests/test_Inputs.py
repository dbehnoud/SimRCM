
import simrcm.sim_tools

#assert inputs.z == 10
def test_bore():
       inputs = simrcm.sim_tools.Inputs('simrcm/tests/example_input.yaml')
       assert inputs.bore == 0.0508 
#assert inputs.v_rcm == 5.47669375000E+002         
#assert inputs.vprofile['vproVol'][2] == 5.43427034574E+002     
        
