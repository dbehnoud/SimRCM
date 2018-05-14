
from .. import Sim
import pytest

def test_sim():

 gas = Sim.simulation()
 
 assert gas.T == 2869.88
 assert gas.P == 259246
