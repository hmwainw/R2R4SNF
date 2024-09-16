###################################################################################
"""
    Last edited on February 19, 2019

    @author: matz

    comments: Test the methods and attributes of the Stage class
    
"""
###################################################################################
import os
import pytest
import numpy as np
import pandas as pd
from collections.abc import Sequence
from nwpy import stage
###################################################################################
# Test fixtures: General cases and special cases
# EG    Fuel cycle                  Stages      Isotopics
#----------------------------------------------------------
# O1    Once-through PWR            1           standard
# 06    Once-through FFH            1           ffh
# 09    Limited-recycle SFR (B&B)   1           batch
# 10    Limited-recycle MSR         1           msr
# 11*   Limited-recycle SFR (B&B)   1           batch
# 13    Limited-recycle PWR-PWR     2           standard
# 15*   Limited-recycle PWR-SFR     2           standard
# 23    Continuous-recycle SFR      1           df/bb
#----------------------------------------------------------
# * indicates case that is not fully implemented, for whatever reason


###################################################################################
# STAGE INSTANTIATION
###################################################################################


def test_init_invalid_stage_number(general):
    """Test that an invalid stage number results in the proper error"""
    
    # actually testing _get_stage_data
    with pytest.raises(ValueError):
        stage.Stage(general.evaluationgroup, 4, test=True) # no eg with 4 stages
    with pytest.raises(ValueError):
        stage.Stage(general.evaluationgroup, 'fake', test=True) # stgn must be num


def test_init_invalid_evaluation_group():
    """Test that an nonexistant evaluation group results in proper error"""
    
    # actually testing _get_stage_data
    with pytest.raises(IOError):
        stage.Stage('fake', 1, test=True) # no eg named 'fake'
    with pytest.raises(IOError):
        stage.Stage('eg41', 1, test=True) # no eg 41


def test_init_paths(general):
    """Test that stage instance file paths all exist"""

    assert os.path.exists(general.wdpath), 'path to working dir failed'
    assert os.path.exists(general.datapath),'path to data dir failed'


###################################################################################
# DISCHARGE STREAMS
###################################################################################


def test_discharge_streams_general(general):
    """Test that the stage can discharge streams for general reactor type"""

    strm = general.discharge_streams()
    assert strm.mass == sum(strm.comp[0.0])


def test_discharge_streams_dfbb(dfbb):
    """Test that the stage can discharge df, bb streams"""

    streams = dfbb.discharge_streams()
    assert len(streams) == 2, 'incorrect number of streams discharged'
    for strm in streams:
        assert strm.mass == sum(strm.comp[0.0])


#def test_discharge_streams_msr(msr):
#    """Test that the stage can discharge msr streams"""
#
#    streams = msr.discharge_streams()
#    assert len(streams) == 2, 'incorrect number of streams discharged'
#    assert hasattr(streams[1], 'index') # list quacks like duck
#    for strm in streams:
#        if(hasattr(strm, 'index')): # list quacks like duck
#            for substrm in strm:
#                assert substrm.mass == sum(substrm.comp.values())
#        else:
#            assert strm.mass == sum(strm.comp.values())


def test_discharge_streams_ffh(ffh):
    """Test that the stage can discharge ffh streams"""

    streams = ffh.discharge_streams()
    assert len(streams) == 2, 'incorrect number of streams discharged'
    for strm in streams:
        assert strm.mass == sum(strm.comp[0.0])


def test_discharge_streams_bnbsfr(bnbsfr):
    """Test that the stage can discharge b&b sfr streams"""

    streams = bnbsfr.discharge_streams()
    for strm in streams:
        if(hasattr(strm, 'index')): # list quacks like duck
            for substrm in strm:
                assert substrm.mass == sum(substrm.comp[0.0])
        else:
            assert strm.mass == sum(strm.comp[0.0])


###################################################################################
# COOLING AND DECAYING STREAMS
###################################################################################


def test_cool_input_nested_lists(general):
    """Test that the cool method can take a strange nested list as input"""

    strm = general.discharge_streams()
    funny_nested_list = [strm, [strm, strm]]
    strm_out = general.cool(funny_nested_list)
    depth = lambda L: isinstance(L, list) and max(map(depth, L))+1
    assert depth(strm_out) == 2
    assert len(strm_out) == 2


def test_cool_times_values(general):
    """Test that the cool method can properly handle different cooling times"""

    strm = general.discharge_streams()
    m0 = strm.mass
    # ZERO CASE
    general.data['coolingtime'] = 0.0 # modify cooling time
    # modify comp to isolate unstable nuclide
    strm.comp = pd.DataFrame({0.0: [strm.mass]}, index=['pu241'])
    strm = general.cool(strm)
    assert 0.0 in strm.comp.columns
    err = (m0-strm.comp[0.0]['pu241'])/m0
    assert err < 0.001, 'zero time cooling gives difference > 0.1%'
    # NEGATIVE CASE
    general.data['coolingtime'] = -1.0 # modify cooling time
    with pytest.raises(AssertionError):
        strm = general.cool(strm)


def test_decay_stream_cool_output_value(general):
    """assert that cooling call to origen returns correct numerical value"""
    
    strm = general.discharge_streams()
    # override the composition
    strm.comp = pd.DataFrame({0.0: [strm.mass]}, index=['pu241'])
    strm.form = 'pu241_test_case'
    # expected mass of pu241 after decay
    ans = strm.mass*np.exp(-np.log(2)*general.data['coolingtime']/14.29)
    strm = general._decay_stream(strm, 'last', 5.0, 5.0, 'I',
                                 cooling=True, rerun=True, remove=True)
    err = (ans - strm.comp[5.0]['pu241'])/ans
    assert err < 0.001, 'difference between origen and numerical decay > 0.1%'


def test_cool_stream_mass_difference(general):
    """assert that cooling call to origen returns correct numerical value"""
    
    strm = general.discharge_streams()
    cooled_strm = general.cool(strm)
    mass_difference = (strm.mass-cooled_strm.mass)/strm.mass
    assert abs(mass_difference) < 0.001, 'mass diff across cooling > 0.1%'


def test_decay_stream_attributes(general):
    """Test that stage.decay method assigns proper attributes"""
    
    strm = general.discharge_streams()
    # expected mass of pu241 after decay
    strm = general.decay(strm)
    assert all([hasattr(strm, a) for a in ['comp', 'heat', 'act', 'tox']])
    

def test_decay_stream_output(general):
    """Test that decay stream gives correct numerical values in attributes"""
    
    time = 40.0
    hl = 14.29
    h = 3.284e-3 # W/g
    strm = general.discharge_streams()
    # override the composition
    strm.comp = pd.DataFrame({0.0: [strm.mass]}, index=['pu241'])
    strm.form = 'pu241_decay_test_case'
    strm = general.decay(strm, endtime=time)
    # expected mass of pu241 after decay
    m = strm.mass*np.exp(-np.log(2)*time/hl)
    n = m*6.02e23/241 # atoms
    lam = np.log(2)/(hl*365*24*3600)
    act = n*lam/3.7e10 # Ci
    heat = h*m
    assert (m - strm.comp[time]['pu241'])/m < 0.001
    assert (act - strm.act[time]['pu241'])/act < 0.001
    assert (heat - strm.heat[time]['pu241'])/heat < 0.001


def test_decay_nonzero_starttime(general):
    """Test that if starttime is nonzero, result times are correct"""

    strm = general.discharge_streams()
    strm = general.cool(strm)
    new_strm = general.decay(strm, starttime='last', endtime=100.0)
    assert max(new_strm.comp.columns)==general.data['coolingtime']+100.0
    newer_strm = general.decay(new_strm, starttime='last', endtime=207.6)
    assert (round(max(newer_strm.comp.columns),2)==
            round(max(new_strm.comp.columns)+207.6,2))


###################################################################################
# REPROCESSING STREAMS
# test the complicated ones
###################################################################################


#def test_reprocess_msr(msr):
#    """Test reprocessing of a complicated, multilevel list of streams"""
#
#    streams = msr.discharge_streams()
#    streams = msr.cool(streams)
#    rw = msr.reprocess(streams) # form: [fp1, [dischargedsalt, fp2]]
#    assert rw[0].mass_fraction('u') == 0.0 # no u in fp from salt treatment
#    assert round(rw[0].mass,2) == round(streams[0].mass,2)
#    for sp in msr.data['wastespecies'].split(','):
#        idx = msr.data['wastespecies'].split(',').index(sp)
#        m = msr.data['wastemass'][idx]
#        assert round(rw[1].mass*rw[1].mass_fraction(sp),2) == m*1e6


def test_reprocess_bnbsfr(bnbsfr):
    """Test reprocessing of a complicated, multilevel list of streams"""

    streams = bnbsfr.discharge_streams()
    streams = bnbsfr.cool(streams)
    u_discharged = 0.0
    for i in range(0, len(streams[0])):
        u_discharged += streams[0][i].mass*streams[0][i].mass_fraction('u')
    rw = bnbsfr.reprocess(streams)
    u_in_waste = 0.0
    for i in range(0, len(rw)-1):
        u_in_waste += rw[i].mass*rw[i].mass_fraction('u')
    assert round(u_discharged*0.01,2) == round(u_in_waste,2) # assume sep eff  99%
    assert streams[len(streams)-1].form == rw[len(rw)-1].form # discharged batch
    assert round(streams[len(streams)-1].mass,2) == round(rw[len(rw)-1].mass,2)


###################################################################################
# LOADING WASTES
###################################################################################
# single stream snf
# multi-stream snf




###################################################################################
# DECAYING WASTES
###################################################################################


###################################################################################
# DISCHARGING WASTES (integration of all of the above)
###################################################################################


# METHODS LEFT TO TEST
# reprocess
# load_waste
# decay
# decay_stream
# write_data
# benchmark_stage
# _sum_stage_waste_streams
# discharge_wastes

