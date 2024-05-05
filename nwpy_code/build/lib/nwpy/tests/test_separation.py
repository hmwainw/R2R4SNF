###############################################################################
"""
    Last edited on February 27, 2019

    @author: matz

    comments: Test the generation and attributes of Separation class
    
"""
###############################################################################
import os
import copy
import pytest
import numpy as np
import pandas as pd
from nwpy import separation
from nwpy import stream
###############################################################################


@pytest.fixture
def echem(dfbb):
    """ """
    print('\n'+'Setup electrochemical solid fuel separation instance')
    return(dfbb.separation)


@pytest.fixture
def aqueous(general_urex):
    """ """
    print('\n'+'Setup aqueous solid fuel separation instance')
    return(general_urex.separation)


@pytest.fixture
def meltrefining(bnbsfr):
    """ """
    print('\n'+'Setup melt-refining solid fuel separation instance')
    return(bnbsfr.separation)


@pytest.fixture
def simple_stream(): # fake composition, fake heat
    print('\n'+'Setup simple Stream instance')
    c = pd.DataFrame({5.0: [1.0, 1.0, 1.0, 5.0, 90.0, 2.0]},
                     index=['sr90','cs137','th232','u235','u238','pu239'])
    h = pd.DataFrame({5.0: [50.0, 20.0, 0.1, 1.0, 0.1, 5.0]},
                     index=['sr90','cs137','th232','u235','u238','pu239'])
    return(stream.Stream(100.0, comp=c, heat=h, id='snf', form='snf',
                         evaluationgroup='test-simple', stagenumber=1))


@pytest.fixture
def msr(msr):
    """ """
    print('\n'+'Setup liquid fuel separation instance')
    return(msr.separation)


@pytest.fixture
def msr_df_stream(): # fake composition, fake heat
    print('\n'+'Setup MSR DF Stream instance')
    c = pd.DataFrame({5.0: [7.1e5, 0.5e5, 2.1e5, 0.3e5]},
                     index=['th232', 'pa233', 'u233', 'u234'])
    h = pd.DataFrame({5.0: [0.1, 4.5, 0.4, 0.2]},
                     index=['th232', 'pa233', 'u233', 'u234'])
    return(stream.Stream(1000000.0, comp=c, id='snf', heat=h, form='df1',
                         evaluationgroup='test-msr-df', stagenumber=1))


#@pytest.fixture
#def msr_fp_stream(): # fake composition, fake heat; normalized to 1e6
#    print('\n'+'Setup MSR FP Stream instance')
#    # gases leave, reactive fp plate out
#    c = pd.DataFrame({0.0: [2.0e5, 1.2e5, 2.8e5, 2.6e5, 1.4e5]},
#                     index=['sr90', 'i129', 'i131', 'cs137', 'ba138'])
#    h = pd.DataFrame({0.0: [10.0, 0.5, 60.1, 4.0, 1.3]},
#                     index=['sr90', 'i129', 'i131', 'cs137', 'ba138'])
#    return(stream.Stream(1000000.0, comp=c, id='snf', heat=h, form='fp2',
#                         evaluationgroup='test-msr-df', stagenumber=1))


###############################################################################
# ELECTROCHEMICAL SEPARATION
###############################################################################


def test_echem_instantiation(echem):
    """assert sep instance has required attributes"""
    
    assert hasattr(echem, 'method')
    assert hasattr(echem, 'recovered')
    assert hasattr(echem, 'data')
    assert hasattr(echem, 'datapath')
    assert hasattr(echem, 'tol')


def test_echem_import_separation_data(echem):
    """test that sep instance can import data"""
    
    file = echem.method+'.sep'
    confirm = open(os.path.join(echem.datapath, 'sep', file)).read()
    confirm = confirm.splitlines()
    test, path = echem._import_separation_data()
    assert os.path.exists(path)
    assert confirm == test


def test_echem_get_output_streams(echem):
    """test that sep instance can read proper wasteforms"""

    confirm = ['CERAMIC', 'GAS', 'METAL'] # echem wasteforms
    data, path = echem._import_separation_data()
    strms, hlines = echem._get_output_streams(data)
    #print wf
    assert all([x in confirm for x in strms])
    assert len(strms) == len(confirm)


def test_echem_read_sep_data(echem):
    """test sep instance can read data and produce dict"""
    
    array, names = echem._read_sep_data()
    assert all([x in ['CERAMIC', 'METAL', 'GAS'] for x in names])
    assert len(array) == 99 # elements up to Es
    assert all(np.sum(array, axis=1) == np.ones(len(array)))


def test_echem_reprocess(echem, simple_stream):
    """test sep instance to apply reprocessing data to stream"""

    # resulting streams should have 0.1% of u, pu (the rest recovered)
    # should have all fp, unrecovered actinides
    outstrm = echem.reprocess(simple_stream, 0.99)
    total = stream.empty()
    for i in range(0, len(outstrm)):
        total = total.add(outstrm[i])
    #outstrm = outstrm[0].add(outstrm[1].add(outstrm[2]))
    assert (round(total.comp[5.0]['sr90'],2) ==
            simple_stream.comp[5.0]['sr90'])
    assert (round(total.comp[5.0]['u235'],2) ==
            0.01*simple_stream.comp[5.0]['u235'])
    assert (round(total.heat[5.0]['pu239'],2) ==
            0.01*simple_stream.heat[5.0]['pu239'])


###############################################################################
# AQUEOUS SEPARATION
###############################################################################


def test_aq_instantiation(aqueous):
    """assert sep instance has required attributes"""
    
    assert hasattr(aqueous, 'method')
    assert hasattr(aqueous, 'recovered')
    assert hasattr(aqueous, 'data')
    assert hasattr(aqueous, 'datapath')
    assert hasattr(aqueous, 'tol')


def test_aq_import_separation_data(aqueous):
    """test that sep instance can import data"""
    
    file = aqueous.method+'.sep'
    confirm = open(os.path.join(aqueous.datapath, 'sep', file)).read()
    confirm = confirm.splitlines()
    test, path = aqueous._import_separation_data()
    assert os.path.exists(path)
    # print(confirm)
    print(test)
    print(dir(test))
    assert confirm == test


def test_aq_reprocess(aqueous, simple_stream):
    """test sep instance to apply reprocessing data to stream"""
    
    # resulting streams should have 5% of u, pu (the rest recovered)
    # should have all fp, unrecovered actinides
    outstrm = aqueous.reprocess(simple_stream, 0.95)
    total = stream.empty()
    for i in range(0, len(outstrm)):
        total = total.add(outstrm[i])
    print(total.comp.columns)
    assert round(total.comp[5.0]['sr90'],2) == simple_stream.comp[5.0]['sr90']
    assert (round(total.comp[5.0]['u235'],2) ==
            0.05*simple_stream.comp[5.0]['u235'])
    assert (round(total.heat[5.0]['pu239'],2) ==
            0.05*simple_stream.heat[5.0]['pu239'])


###############################################################################
# MELT REFINING
###############################################################################


def test_mr_recovery_dict(meltrefining):
    """test that melt refining sep makes special recovery dict"""

    rd = meltrefining._build_recovery_dict(0.99)
    print(meltrefining.recovered)
    assert len(rd) > len(meltrefining.recovered.split(','))
    assert 'am' in rd.keys()
    assert rd['am']==0.05 #hardcoded


def test_mr_custom_rdict(meltrefining):
    """Test that user can provide custom separations criteria"""

    rd=meltrefining._build_recovery_dict(0.99,user_sp_dict={'pu':0.8,'am':0.7})
    assert all([x in rd.keys() for x in ['pu', 'am']])
    assert rd['pu']==0.8
    assert rd['am']==0.7


###############################################################################
# MSR SEPARATIONS
###############################################################################


def test_msr_recovery_hilo(msr):
    """test that msr sep makes special recovery dict"""
    
    rd_lo = msr._build_recovery_dict(0.99, re='lo')
    rd_hi = msr._build_recovery_dict(0.99, re='hi')
    assert round(1-rd_lo['eu'], 3) == round((1-rd_hi['eu'])/10.0,3)
    assert round(1-rd_lo['gd'], 3) == round((1-rd_hi['gd'])/10.0,3)
    assert round(1-rd_lo['sm'], 3) == round((1-rd_hi['sm'])/10.0,3)


def test_msr_custom_rdict(msr, msr_df_stream):
    """Test that user can provide custom separations criteria"""

    input = {'th':0.915}
    rd = msr._build_recovery_dict(0.97, user_sp_dict=input)
    assert rd['th']==input['th']
    w = msr.reprocess(msr_df_stream, 0.97, input)
    u_in = msr_df_stream.mass*msr_df_stream.mass_fraction('u')
    th_in = msr_df_stream.mass*msr_df_stream.mass_fraction('th')
    assert round(w.mass_fraction('u')*w.mass,3) == round((1-0.97)*u_in,3)
    assert round(w.mass_fraction('th')*w.mass,3) == round((1-input['th'])*th_in,3)


###############################################################################

