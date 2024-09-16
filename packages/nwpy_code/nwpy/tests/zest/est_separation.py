###############################################################################
"""
    Last edited on September 19, 2018

    @author: matz

    comments: Test the generation and attributes of Separation class
    
"""
###############################################################################
import os
import copy
import pytest
import numpy as np
import pandas as pd
from nwpy.fuelcycle import separation
from nwpy.fuelcycle import stream
###############################################################################
# FIXTURES: SolidFuelSep, LiquidFuelSep
###############################################################################


@pytest.fixture
def solid_fuel_sep(dfbb):
    """ """
    print('\n'+'Setup solid fuel separation instance')
    return(dfbb.separation)


@pytest.fixture
def simple_stream(): # fake composition, fake heat
    print('\n'+'Setup simple Stream instance')
    c = {'sr90':1.0, 'cs137':1.0, 'th232':1.0,
         'u235':5.0, 'u238':90.0, 'pu239':2.0}
    h = {'sr90':50.0, 'cs137':20.0, 'th232':0.1,
        'u235':1.0, 'u238':0.1, 'pu239':5.0}
    return(stream.Stream(100.0, c, id='snf', heat=h, form='snf', time=5.0,
                         evaluationgroup='test-simple', stagenumber=1))


@pytest.fixture
def liquid_fuel_sep(msr):
    """ """
    print('\n'+'Setup liquid fuel separation instance')
    return(msr.separation)


@pytest.fixture
def msr_df_stream(): # fake composition, fake heat
    print('\n'+'Setup MSR DF Stream instance')
    c = {'th232':7.1e5, 'pa233':0.5e5, 'u233':2.1e5, 'u234':0.3e5}
    h = {'th232':0.1, 'pa233':4.5, 'u233':0.4, 'u234':0.2}
    return(stream.Stream(sum(c.values()), c, id='snf', heat=h, form='df1',
                         evaluationgroup='test-msr-df', stagenumber=1,
                         time=5.0))


@pytest.fixture
def msr_fp_stream(): # fake composition, fake heat; normalized to 1e6
    print('\n'+'Setup MSR FP Stream instance')
    # gases leave, reactive fp plate out
    c = {'sr90':2.0e5,'i129':1.2e5,'i131':2.8e5,'cs137':2.6e5,'ba138':1.4e5}
    h = {'sr90':10.0, 'i129':0.5, 'i131':60.1, 'cs137':4.0, 'ba138':1.3}
    return(stream.Stream(sum(c.values()), c, id='hlw', heat=h, form='fp2',
                         evaluationgroup='test-msr-fp', stagenumber=1,
                         time=5.0))


###############################################################################
# SOLID FUEL SEPARATION (sfs)
###############################################################################


def test_sfs_instantiation(solid_fuel_sep):
    """assert sep instance has required attributes"""
    
    assert hasattr(solid_fuel_sep, 'method')
    assert hasattr(solid_fuel_sep, 'recovered')
    assert hasattr(solid_fuel_sep, 'data')
    assert hasattr(solid_fuel_sep, 'datapath')
    assert hasattr(solid_fuel_sep, 'tol')


def test_sfs_import_separation_data(solid_fuel_sep):
    """test that sep instance can import data"""
    
    file = (solid_fuel_sep.method+'_'+
            '-'.join(solid_fuel_sep.recovered.split(','))+'_fces.sep')
    confirm = open(os.path.join(solid_fuel_sep.datapath, 'sep/'+file)).read()
    confirm = confirm.splitlines()
    test, path = solid_fuel_sep._import_separation_data()
    assert os.path.exists(path)
    assert confirm == test


def test_sfs_get_output_streams(solid_fuel_sep):
    """test that sep instance can read proper wasteforms"""

    confirm = ['CERAMIC', 'GAS', 'METAL', 'U', 'PU'] # echem wasteforms
    data, path = solid_fuel_sep._import_separation_data()
    strms, ids, hlines = solid_fuel_sep._get_output_streams(data)
    #print wf
    assert all([x in confirm for x in strms])
    assert len(strms) == len(confirm)


def test_sfs_read_sep_data(solid_fuel_sep):
    """test sep instance can read data and produce dict"""
    
    testdict = solid_fuel_sep._read_sep_data(False)
    assert all([x in ['CERAMIC', 'METAL', 'GAS'] for x in testdict['streams']])
    assert len(testdict['eff']) == 99 # elements up to Es
    # try with products on
    testdict = solid_fuel_sep._read_sep_data(True)
    assert all([x in ['CERAMIC', 'METAL', 'GAS', 'U', 'PU']
                for x in testdict['streams']])


def test_sfs_reprocess(solid_fuel_sep, simple_stream):
    """test sep instance to apply reprocessing data to stream"""

    # resulting streams should have 0.1% of u, pu (the rest recovered)
    # should have all fp, unrecovered actinides
    outstrm = solid_fuel_sep.reprocess(simple_stream)
    outstrm = outstrm[0]+outstrm[1]+outstrm[2]
    assert outstrm.comp['sr90'] == simple_stream.comp['sr90']
    assert outstrm.comp['u235'] == 0.01*simple_stream.comp['u235']
    assert outstrm.heat['pu239'] == 0.01*simple_stream.heat['pu239']


###############################################################################
# LIQUID FUEL SEPARATION (lfs)
###############################################################################


def test_lfs_instantiation(liquid_fuel_sep):
    """assert sep instance has required attributes"""
    
    assert hasattr(liquid_fuel_sep, 'wastespecies')
    assert hasattr(liquid_fuel_sep, 'dischargedfuel')
    assert hasattr(liquid_fuel_sep, 'wastemass')
    assert hasattr(liquid_fuel_sep, 'data')
    assert hasattr(liquid_fuel_sep, 'datapath')
    assert hasattr(liquid_fuel_sep, 'tol')


def test_lfs_apply_tol(liquid_fuel_sep, simple_stream):
    """test sep instance application of tolerance"""

    liquid_fuel_sep.tol = 2.0 # silly
    temp_stream = copy.deepcopy(simple_stream)
    new_stream = liquid_fuel_sep._apply_tol(temp_stream)
    assert not any([x in  new_stream.comp for x in ['sr90', 'cs137', 'th232']])
    print simple_stream.comp.keys()
    loss = (simple_stream.comp['sr90']+simple_stream.comp['cs137']+
            simple_stream.comp['th232'])
    assert new_stream.mass == simple_stream.mass - loss


def test_lfs_get_fp_stream_mass(liquid_fuel_sep, msr_fp_stream):
    """test the get_fp_stream_mass function"""

    temp_stream = copy.deepcopy(msr_fp_stream)
    new_stream = liquid_fuel_sep._get_fp_stream_mass(temp_stream)
    mass_fp=liquid_fuel_sep.wastemass[liquid_fuel_sep.wastespecies.index('fp')]
    assert round(new_stream.mass,2) == round(mass_fp*1e6,2)


def test_lfs_makeup(liquid_fuel_sep, msr_fp_stream, msr_df_stream):
    """test the makeup function"""
    
    temp_fp = copy.deepcopy(msr_fp_stream)
    temp_df = copy.deepcopy(msr_df_stream)
    temp_fp = liquid_fuel_sep._get_fp_stream_mass(temp_fp)
    temp_fp = liquid_fuel_sep._makeup(temp_fp, temp_df)
    mass_th=liquid_fuel_sep.wastemass[liquid_fuel_sep.wastespecies.index('th')]
    assert round(temp_fp.comp['th232'],2) == mass_th*1e6
    mass_u=liquid_fuel_sep.wastemass[liquid_fuel_sep.wastespecies.index('u')]
    assert round(temp_fp.comp['u233']+temp_fp.comp['u234'],2) == mass_u*1e6

            
def test_lfs_reprocess(liquid_fuel_sep, msr_fp_stream, msr_df_stream):
    """test the reprocess fxn"""

    outstreams = liquid_fuel_sep.reprocess([msr_df_stream, msr_fp_stream])
    for strm in outstreams:
        if(strm.form =='fp2'):
            for i in xrange(0, len(liquid_fuel_sep.wastespecies)):
                ws = liquid_fuel_sep.wastespecies[i]
                m = liquid_fuel_sep.wastemass[i]*1e6
                assert round(strm.mass_fraction(ws)*strm.mass,2) == m


###############################################################################