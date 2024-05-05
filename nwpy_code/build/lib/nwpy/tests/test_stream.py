###############################################################################
"""
    Last edited on February 19, 2019

    @author: matz

    comments: Test the generation and operations of the Stream instance

"""
###############################################################################
import pytest
import numpy as np
import pandas as pd
from copy import deepcopy
from nwpy import stream
###############################################################################
# STILL TO DO
# Stream subtraction


def test_stream_init_with_invalid_mass():
    """Test that initialization w incorrect mass returns error"""
    data = {0.0: 0.0}
    with pytest.raises(AssertionError): # negative mass
        stream.Stream(-1.0,comp=pd.DataFrame(data,index=['a']),form='test')
    with pytest.raises(AssertionError):  # non-number mass
        stream.Stream('fake',comp=pd.DataFrame(data,index=['a']),form='test')


def test_strm_addition():
    """Test the combination of two streams"""
    
    m1, m2 = np.random.random(2)*100.0
    c1 = pd.DataFrame({0.0: [np.random.random()*10.0,
                             np.random.random()*10.0]},
                      index=['a', 'b'])
    c2 = pd.DataFrame({0.0: [np.random.random()*10.0,
                             np.random.random()*10.0]},
                      index=['a', 'c'])
    h1 = pd.DataFrame({0.0: [np.random.random(),
                             np.random.random()]},
                      index=['a', 'b'])
    h2 = pd.DataFrame({0.0: [np.random.random(),
                             np.random.random()]},
                      index=['a', 'c'])
    s1 = stream.Stream(m1, comp=c1, heat=h1, form='test1')
    s2 = stream.Stream(m2, comp=c2, heat=h2, form='test2')
    s_test = s1.add(s2)
    assert s_test.mass == m1+m2
    assert s_test.comp[0.0]['a'] == c1[0.0]['a']+c2[0.0]['a']
    assert s_test.comp[0.0]['b'] == c1[0.0]['b']
    assert s_test.comp[0.0]['c'] == c2[0.0]['c']
    assert s_test.heat[0.0]['a'] == h1[0.0]['a']+h2[0.0]['a']
    assert s_test.heat[0.0]['b'] == h1[0.0]['b']
    assert s_test.heat[0.0]['c'] == h2[0.0]['c']
    
    
def test_strm_add_mismatched_times():
    """Test that assertion error catches streams with mismatched times"""
    
    m1, m2 = np.random.random(2)*100.0
    c1 = pd.DataFrame({0.0: [np.random.random()*10.0]}, index=['a'])
    c2 = pd.DataFrame({1.0: [np.random.random()*10.0]}, index=['a'])
    h1 = pd.DataFrame({0.0: [np.random.random()]}, index=['a'])
    h2 = pd.DataFrame({1.0: [np.random.random()]}, index=['a'])
    s1 = stream.Stream(m1, comp=c1, heat=h1, form='test1')
    s2 = stream.Stream(m2, comp=c2, heat=h2, form='test2')
    with pytest.raises(AssertionError):
        s_test = s1.add(s2, time='last')
    with pytest.raises(AssertionError):
        s_test = s1.add(s2, time=0.0)


def test_strm_add_empty_strm():
    m1=np.random.random()*100.0
    c1=pd.DataFrame({0.0:[np.random.random()*10.0,np.random.random()*10.0]},
                    index=['a', 'b'])
    h1 = pd.DataFrame({0.0:[np.random.random(),np.random.random()]},
                      index=['a', 'b'])
    s1 = stream.Stream(m1, comp=c1, heat=h1, form='test1')
    assert s1 == stream.empty().add(s1)
    assert s1 == s1.add(stream.empty())


def test_mass_fraction_inputs():
    """Test that the mass fraction method works for different inputs"""
    
    # Isotopes
    c = {0.0: [5.0, 95.0]}
    idx = ['u235', 'u238']
    s = stream.Stream(100.0,comp=pd.DataFrame(c,index=idx),form='test')
    assert s.mass_fraction('u235') == 0.05
    # Elements
    c = {0.0: [5.0, 93.0, 2.0]}
    idx = ['u235', 'u238', 'pu239']
    s = stream.Stream(100.0,comp=pd.DataFrame(c,index=idx),form='test')
    assert s.mass_fraction('u') == 0.98
    assert s.mass_fraction('pu') == 0.02
    # Groups/species
    c = {0.0: [1.5, 1.0, 2.0, 60.0, 20.0, 5.0, 4.0, 6.5]}
    idx = ['sr90','cs134','ba137','th232','u238','np237','pu239','am241m']
    s = stream.Stream(100.0,comp=pd.DataFrame(c,index=idx),form='test')
    assert s.mass_fraction('act') == 0.955
    assert s.mass_fraction('fp') == 0.045
    assert s.mass_fraction('ma') == 0.115
    assert s.mass_fraction('tru') == 0.155


###############################################################################
