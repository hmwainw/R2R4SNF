###############################################################################
"""
    Last edited on February 19, 2019

    @author: matz

    comments: Tests for the heat module

    
"""
###############################################################################
import os
import math
import pytest
import numpy as np
from copy import deepcopy
from nwpy import stream
from nwpy.repository_area import heat
###############################################################################
# integrate
# _point
# _finiteline
# wall temp
# wall temp single source
# wp temp
#

###############################################################################


@pytest.fixture
def ht():
    print('\n'+'Setup heat transfer instance; surface storage time = 0.0')
    return(heat.HeatTransfer())


###############################################################################


def test_ht_instantiation(ht):
    """Test that ht instance has proper attributes"""

    assert ht.st == 0.0
    assert ht.tol == 1.0e-4


def test_ht_integration_simple(ht):
    
    t1 = int(np.random.random()*(6-3)+3) # random int between 3 and 6
    c = np.random.random()*(7-2)+2 # random number between 2 and 7
    ans = round((t1**2)/2.+c*t1,3)
    def f(x,c):
        return(x+c)
    assert round(ht._integrate(f, t1, arguments=(c,)),3)==ans


def test_point_src_constant_heat(granite_const_heat, ht):
    """Test that integration with constant heat generation gives 
    correct analytical result"""
        
    r = deepcopy(granite_const_heat)
    # set up other stuff
    trgt = float(int(int(np.random.random()*(100-10)+10))) # random time (y)
    dist = int(np.random.random()*(20-2)+2) # random distance [2,20] (m)
    rad = math.sqrt(dist**2+r.ebs['r_drift']**2)
    a = r.rock['a']*365*24*3600
    k = r.rock['k']
    # const heat, eqn reduces to erfc
    ans = 10.0/4.0/np.pi/k/rad
    ans *= math.erfc(rad/math.sqrt(4.0*a*trgt))
    assert round(ht._point(trgt, r, dist),5)==round(ans,5)


def test_finite_line_src_constant_heat(granite_const_heat, ht):
    """Test that integration with constant heat generation gives 
    correct result"""
        
    r = deepcopy(granite_const_heat)
    trgt = float(int(int(np.random.random()*(100-10)+10))) # y
    dr = r.ebs['r_drift']
    a = r.rock['a']*365*24*3600
    k = r.rock['k']
    wp_len = r.pkg['l']
    # const heat
    c = 10.0/8.0/np.pi/k/wp_len
    # still need to integrate b/c of erf fxns
    def integrand(t):
        denom = trgt-t
        # no heat term
        expterm = np.exp(-(dr**2)/4.0/a/(trgt-t))
        erf1 = math.erf(0.25*wp_len/np.sqrt(a*(trgt-t)))
        erf2 = math.erf(-0.25*wp_len/np.sqrt(a*(trgt-t)))
        return(expterm*(erf1-erf2)/denom)
    # "semi-analytical" ans
    ans = c*ht._integrate(integrand, trgt)
    assert round(ht._finiteline(trgt, r),5)==round(ans,5)


def test_calc_wall_temp_const_heat(granite_const_heat, ht):
    """Test that overlaying many sources with constant heat 
    gives the correct result"""
        
    r = deepcopy(granite_const_heat)
    trgt = float(int(int(np.random.random()*(100-10)+10))) # y
    Twall = 0.0
    # decrease array size so just nine total sources
    r.N = 3 # on each side; 3x3 total
    # other stuff
    dr = dr = r.ebs['r_drift']
    a = r.rock['a']*365*24*3600
    k = r.rock['k']
    wp_len = r.pkg['l']
    # CENTRAL PACKAGE
    # still need to integrate b/c of erf fxns
    def integrand(t):
        denom = trgt-t
        # no heat term
        expterm = np.exp(-(dr**2)/4.0/a/(trgt-t))
        erf1 = math.erf(0.25*wp_len/np.sqrt(a*(trgt-t)))
        erf2 = math.erf(-0.25*wp_len/np.sqrt(a*(trgt-t)))
        return(expterm*(erf1-erf2)/denom)
    Twall += (10.0/8.0/np.pi/k/wp_len)*ht._integrate(integrand, trgt)
    # adjacent packages in central drift (2)
    d = math.sqrt(r.spacing['pkg']**2+r.ebs['r_drift']**2)
    Twall += 2*(10.0/4.0/np.pi/k/d)*math.erfc(d/math.sqrt(4.0*a*trgt))
    # packages in adjacent drift (3); two corner packages (in both drifts)
    d=math.sqrt(r.spacing['pkg']**2+r.spacing['drift']**2+r.ebs['r_drift']**2)
    Twall += 4*(10.0/4.0/np.pi/k/d)*math.erfc(d/math.sqrt(4.0*a*trgt))
    # packages in adjacent drift (3); two packages directly adjacent
    d = math.sqrt(r.spacing['drift']**2+r.ebs['r_drift']**2)
    Twall += 2*(10.0/4.0/np.pi/k/d)*math.erfc(d/math.sqrt(4.0*a*trgt))
    ans = Twall+r.ambient_temp
    assert round(ht._calc_wall_temp(trgt, r),5) == round(ans,5)


def test_calc_wall_temp_single_source_const_heat(granite_const_heat, ht):
    """Test the result returned for only finite line source"""
    
    r = deepcopy(granite_const_heat)
    trgt = float(int(np.random.random()*(100-10)+10)) # y
    ans = ht._finiteline(trgt, r)+r.ambient_temp
    assert round(ht._calc_wall_temp(trgt,r,allsources=False),5) == round(ans,5)


def test_calc_wp_temp_const_heat(granite_const_heat, ht):
    """Test that inner model calculation gives correct 
    result for simple case"""
    
    r = deepcopy(granite_const_heat)
    trgt = float(int(np.random.random()*(100-10)+10)) # y
    # overwrite all thermal conductivities with "1"
    r.ebs['k'] = np.ones(len(r.ebs['k'])-1)
    # assume some number for the outer temperature
    T_out = float(int(np.random.random()*(100-20)+20))
    # generate test answer
    ans = T_out
    print(ans)
    r_out = r.ebs['r_drift']
    for l in range(0, len(r.ebs['k'])):
        r_in = r_out-r.ebs['dr'][l]
        ans += r.decay_heat(1)*np.log(r_out/r_in)/2.0/r.pkg['l']/np.pi
    assert round(ht._calc_wp_temp(trgt, r, T_out), 5)==round(ans,5)


###############################################################################
