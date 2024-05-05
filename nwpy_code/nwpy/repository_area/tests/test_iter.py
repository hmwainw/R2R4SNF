###############################################################################
"""
    Last edited on February 19, 2019

    @author: matz

    comments: Tests for the Iter class
    
"""
###############################################################################
import os
import pytest
import numpy as np
from copy import deepcopy
from nwpy.repository_area import iter
###############################################################################


@pytest.fixture
def new_iter():
    print('\n'+'Setup new iter.Iteration instance')
    return(iter.Iteration())


def test_iter_instantiation(new_iter):
    """Test that the new_iter attributes are instantiated correctly"""
    
    assert new_iter.array_idx == 0.0
    assert new_iter.iter_idx == 0.0
    assert len(new_iter.data) == 4
    assert all([v == [] for v in new_iter.data.values()])


def test_iter_update(new_iter):
    """Test that updating the iteration index properly changes attributes"""

    temp = deepcopy(new_iter)
    temp.update()
    assert temp.iter_idx == new_iter.iter_idx+1
    assert temp.array_idx == new_iter.array_idx


def test_iter_read(new_iter):
    """Test that the iter class can read and update internal data"""

    temp = deepcopy(new_iter)
    x = {'wps': 10.0, 'ds': 10.0, 'area': 120.0, 'temp': 130.5}
    temp.read(wps=x['wps'], ds=x['ds'], area=x['area'], temp=x['temp'])
    assert all([temp.data[k][0]==x[k] for k in x.keys()])
    with pytest.raises(KeyError):
        temp.read(fake=5.0)


def test_iter_reset(new_iter):
    """Test that the iter class can reset after reading data"""
    
    temp = deepcopy(new_iter)
    i=0
    while(i<5): # read in some fake data, just gonna delete it
        v = np.random.random(3)
        x = {'wps': v[0]+9.0, 'ds': v[1]+9.0,
            'area': (v[0]+9.0)*(v[1]+9.0), 'temp': v[2]*5+98.0}
        temp.read(wps=x['wps'], ds=x['ds'], area=x['area'], temp=x['temp'])
        i+=1
    temp.reset()
    assert all([v == [] for v in new_iter.data.values()])
    assert temp.array_idx == new_iter.array_idx+1
    assert temp.iter_idx == 0.0


def test_iter_get_peak_temp(new_iter):
    """Test that the iter class can match temperature value with target area"""
    
    temp = deepcopy(new_iter)
    i=0
    while(i<5): # read in some fake data
        v = np.random.random(3)
        x = {'wps': v[0]+9.0, 'ds': v[1]+9.0,
            'area': (v[0]+9.0)*(v[1]+9.0), 'temp': v[2]*5+98.0}
        temp.read(wps=x['wps'], ds=x['ds'], area=x['area'], temp=x['temp'])
        i+=1
    target_area = min(temp.data['area'])
    target_temp = temp.get_peak_temp(target_area)
    test_idx = list(temp.data['area']).index(target_area)
    assert target_temp==temp.data['temp'][test_idx]


###############################################################################
