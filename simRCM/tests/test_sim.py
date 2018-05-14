
from ..simrcm import sim
import pytest

def test_sim():

 gas = Sim.simulation()
 
 assert gas.T == 2869.88
 assert gas.P == 259246
