###############################################################################
"""
    Last edited on February 19, 2019

    @author: matz

    comments: Test the generation and attributes of Reactor class
    
"""
###############################################################################
import os
import pytest
import numpy as np
import pandas as pd
from nwpy import reactor
from nwpy import nuclidedata
###############################################################################
# Test fixtures: General cases and special cases
# These are imported from the Stage fixtures defined in conftest.py
#------------------------------------------------------------------------------


@pytest.fixture
def rx_general(general):
    """Discharge col header 'ST#-discharged (g)' """
    
    print('\n'+'Setup general reactor instance')
    return(general.reactor)


@pytest.fixture
def rx_dfbb(dfbb):
    """Discharge col header 'ST#-('BB' or 'DF')-discharged (g)' """

    print('\n'+'Setup DF-BB reactor instance')
    return(dfbb.reactor)


@pytest.fixture
def rx_ffh(ffh):
    """FFH reactor, special col headers"""

    print('\n'+'Setup FFH reactor instance')
    return(ffh.reactor)


@pytest.fixture
def rx_msr(msr):
    """MSR reactor, special col headers"""
    
    print('\n'+'Setup MSR reactor instance')
    return(msr.reactor)


@pytest.fixture
def rx_bnbsfr(bnbsfr):
    """BnB SFR reactor, special col headers"""
    
    print('\n'+'Setup BnB SFR reactor instance')
    return(bnbsfr.reactor)


###############################################################################
# GENERAL CASE
###############################################################################


def test_general_import_isotopic_csv(rx_general):
    """Test that reactor instance can import csv with isotopic data"""
    
    df = rx_general._import_isotopic_csv()
    file = rx_general.datapath+'iso/'+rx_general.evaluationgroup+'.csv'
    # generate new file by alternate method to confirm
    confirm = open(file).read().splitlines()
    # test the size of the imported array
    assert len(confirm)-1==len(df)
    assert len(confirm[0].split(','))-1==len(df.columns)


def test_general_get_time(rx_general):
    """Test the discharge time stamp function"""
    
    t = rx_general.get_time()
    assert t == 0.0
    # change time, make sure result is still right
    rx_general.data['coolingtime']=0.0
    t = rx_general.get_time()
    assert t == 5.0


def test_general_get_composition_data(rx_general):
    """Test the get_composition_data function"""
    
    df = rx_general._import_isotopic_csv()
    cols, masses, labels = rx_general._streams(df=df)
    comp = rx_general.get_composition_data()
    assert len(comp) == len(cols)
    assert all([x in comp.keys() for x in labels])
    # will fail if bnbsfr, "_df" appended to label
    
    
def test_general_discharge_streams_general(rx_general):
    """Test that each produced stream has values that sum to the total mass"""
    
    outstreams = rx_general.discharge_streams()
    if(not hasattr(outstreams, 'index')):
        outstreams = [outstreams]
    for strm in outstreams:
        assert round(strm.mass, 2) == round(sum(strm.comp[0.0]))


def test_general_make_outstreams(rx_general):
    """Test that make_outstreams properly handles lists and non-lists"""
    
    l = ['a', 'b', 'c', 'd']
    assert hasattr(rx_general._make_outstreams(l), 'index')
    l = ('a', 'b', 'c', 'd')
    assert hasattr(rx_general._make_outstreams(l), 'index')
    l = 'stream'
    assert isinstance(rx_general._make_outstreams(l), str)


def test_general_instantiation(rx_general):
    """Test instantiation and attribute inventory of general reactor class"""
    
    assert hasattr(rx_general, 'data')
    assert hasattr(rx_general, 'datapath')
    assert hasattr(rx_general, 'evaluationgroup')
    assert hasattr(rx_general, 'number')
    assert os.path.exists(rx_general.datapath)


def test_general_streams(rx_general):
    """Test that the _streams method correctly reads df"""
    
    df = rx_general._import_isotopic_csv()
    tup = rx_general._streams(df=df)
    # generate new file by alternate method to confirm
    file = rx_general.datapath+'iso/'+rx_general.evaluationgroup+'.csv'
    confirm = open(file).read().splitlines()
    cols = confirm[0]
    cols = [c for c in cols.split(',') if 'discharge' in c and '1' in c]
    assert len(tup)==3
    assert tup[0] == cols
    assert tup[1] == [rx_general.data['mass']]
    assert tup[2] == ['snf']


###############################################################################
# DF-BB CASE
###############################################################################


def test_dfbb_instantiation(rx_dfbb):
    """Test instantiation and attribute inventory of general reactor class"""
    
    assert hasattr(rx_dfbb, 'data')
    assert hasattr(rx_dfbb, 'datapath')
    assert hasattr(rx_dfbb, 'evaluationgroup')
    assert hasattr(rx_dfbb, 'number')
    assert os.path.exists(rx_dfbb.datapath)


def test_dfbb_streams(rx_dfbb):
    """Test that the _streams method correctly reads df"""
    
    df = rx_dfbb._import_isotopic_csv()
    tup = rx_dfbb._streams(df=df)
    # generate new file by alternate method to confirm
    file = rx_dfbb.datapath+'iso/'+rx_dfbb.evaluationgroup+'.csv'
    confirm = open(file).read().splitlines()
    cols = confirm[0]
    cols = [c for c in cols.split(',') if 'discharge' in c and '1' in c]
    assert len(tup)==3
    assert tup[0] == cols
    assert tup[1] == [rx_dfbb.data['driver'], rx_dfbb.data['blanket']]
    assert tup[2] == ['driver', 'blanket']


###############################################################################
# FFH CASE
###############################################################################


def test_ffh_instantiation(rx_ffh):
    """Test instantiation and attribute inventory of FFH reactor class"""
    
    assert hasattr(rx_ffh, 'data')
    assert hasattr(rx_ffh, 'datapath')
    assert hasattr(rx_ffh, 'evaluationgroup')
    assert hasattr(rx_ffh, 'number')
    assert os.path.exists(rx_ffh.datapath)


def test_ffh_streams(rx_ffh):
    """Test that the _streams method correctly reads df"""
    
    df = rx_ffh._import_isotopic_csv()
    tup = rx_ffh._streams(df=df)
    # generate new file by alternate method to confirm
    file = rx_ffh.datapath+'iso/'+rx_ffh.evaluationgroup+'.csv'
    confirm = open(file).read().splitlines()
    cols = confirm[0].split(',')[1:]
    #cols = [c for c in cols.split(',') if 'discharge' in c and '1' in c]
    assert len(tup)==3
    assert tup[0] == cols
    assert tup[1] == [rx_ffh.data['dischargedfuel'], rx_ffh.data['salttreatment']]
    assert tup[2] == ['snf', 'fp1']


###############################################################################
# MSR CASE
###############################################################################


def test_msr_instantiation(rx_msr):
    """Test instantiation and attribute inventory of MSR reactor class"""
    
    assert hasattr(rx_msr, 'data')
    assert hasattr(rx_msr, 'datapath')
    assert hasattr(rx_msr, 'evaluationgroup')
    assert hasattr(rx_msr, 'number')
    assert os.path.exists(rx_msr.datapath)


def test_msr_streams(rx_msr):
    """Test that the msr _streams method correctly reads df"""
    
    df = rx_msr._import_isotopic_csv()
    tup = rx_msr._streams(df=df)
    file = rx_msr.datapath+'iso/'+rx_msr.evaluationgroup+'.csv'
    confirm = open(file).read().splitlines()
    cols = confirm[0].split(',')[1:]
    assert len(tup)==3
    assert all(len(x)==2 for x in tup) # all should be len 2: fuelsalt, fp1


###############################################################################
# BNB SFR CASE
###############################################################################


def test_bnbsfr_instantiation(rx_bnbsfr):
    """Test instantiation and attribute inventory of BnB SFR reactor class"""
    
    assert hasattr(rx_bnbsfr, 'data')
    assert hasattr(rx_bnbsfr, 'datapath')
    assert hasattr(rx_bnbsfr, 'evaluationgroup')
    assert hasattr(rx_bnbsfr, 'number')
    assert os.path.exists(rx_bnbsfr.datapath)


def test_bnbsfr_streams(rx_bnbsfr):
    """Test that the bnb sfr _streams method correctly reads df"""
    
    df = rx_bnbsfr._import_isotopic_csv()
    tup = rx_bnbsfr._streams(df=df)
    file = rx_bnbsfr.datapath+'iso/'+rx_bnbsfr.evaluationgroup+'.csv'
    confirm = open(file).read().splitlines()
    cols = confirm[0].split(',')[1:]
    #cols = [c for c in cols.split(',') if 'discharge' in c and '1' in c]
    assert len(tup)==3
    assert all(x in cols for x in tup[0])


#def test_bnbsfr_make_outstreams(bnbsfr):
#    raise NotImplementedError


###############################################################################
