###############################################################################
"""
    Last edited on May 6, 2019

    @author: matz

    comments: Tests for the repository module of the repository_area package

    
"""
###############################################################################
import os
from nwpy.repository_area import repository
###############################################################################


def test_instantiation_without_waste(granite_no_waste):
    """Test that all attributes are available and that 
    instantiation works without waste specification"""
    
    assert hasattr(granite_no_waste, 'st')
    assert hasattr(granite_no_waste, 'rock')
    #assert hasattr(granite_no_waste, 'wdpath')
    assert hasattr(granite_no_waste, 'ht')
    assert hasattr(granite_no_waste, 'iter')
    assert hasattr(granite_no_waste, 'datapath')


def test_instantiation_with_waste(granite_const_heat):
    """Test that waste loading is okay when 
    waste is specified at instantiation"""

    assert hasattr(granite_const_heat, 'waste')
    assert hasattr(granite_const_heat, 'decay_heat')
    assert hasattr(granite_const_heat, 'pkg')
    assert all(x in granite_const_heat.waste.keys() for x in ['n_wf', 'id'])
    assert granite_const_heat.decay_heat(500.0) == 10.0


###############################################################################
