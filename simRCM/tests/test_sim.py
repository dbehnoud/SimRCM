from .. import reaction_chamber
import pytest

def test_sim():
    
    """
    check return for two datapoints
    """

    gas = reaction_chamber.simulation()
     
    assert gas.T == pytest.approx(2869.88, rel=1e-1)
    assert gas.P == pytest.approx(259246, rel=1e-1)
