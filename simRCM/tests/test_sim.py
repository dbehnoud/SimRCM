from .. import sim
import pytest

def test_sim():

 gas = sim.simulation()
 
 assert gas.T == pytest.approx(2869.88, rel=1e-1)
 assert gas.P == pytest.approx(259246, rel=1e-1)
