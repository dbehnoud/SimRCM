import simrcm.simulation
import simrcm.sim_tools
import pytest


def test_1():
    ig1, p1, temp1, t1 = simrcm.simulation.simulation1('simrcm/tests/example_input.yaml')
    assert ig1 == pytest.approx(0.03, rel=1e-1)
    
