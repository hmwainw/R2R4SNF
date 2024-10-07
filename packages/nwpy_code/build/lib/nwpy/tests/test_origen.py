###############################################################################
"""
    Last edited on May 22, 2019

    @author: matz

    comments: Test the generation and attributes of Origen class
    
"""
###############################################################################
import os
import numpy as np
import pandas as pd
import subprocess
import pytest
from nwpy import origen
from nwpy import stream
###############################################################################
# Not yet unit tested:
# Origen.update_stream
# Origen.run_origen (sorta tested with test_scale_path)
# Origen.make_file
# Origen.write_temp
###############################################################################


@pytest.fixture
def empty_stream():
    print('\n'+'Setup empty Stream instance')
    return(stream.empty())


@pytest.fixture
def simple_stream():
    print('\n'+'Setup simple Stream instance')
    c = pd.DataFrame({0.0: [1.0]}, index=['pu241'])
    return(stream.Stream(1.0, comp=c, id='snf', form='snf'))


@pytest.fixture
def origen_inst():
    print('\n'+'Setup Origen instance')
    return(origen.Origen)


###############################################################################


def test_update_stream(general, empty_stream, origen_inst):
    """Test the update stream method"""
    
    # write test plt file to working dir
    name = 'test_opus_output.plt'
    infile = open(general.datapath+name).read()
    testfile = open(general.wdpath+'/'+name, 'w')
    testfile.write(infile)
    testfile.close()
    # update empty stream
    strm = origen_inst.update_stream(general.wdpath,empty_stream,name,0.0)
    #assert round(strm.mass, 2) == 100.0
    assert round(sum(strm.comp[5.0]), 2) == 100.0
    assert round(strm.comp[5.0]['u235'],2) == 2.0
    assert round(strm.comp[5.0]['pu241'],2) == 0.61


def test_make_write_file(general, simple_stream, origen_inst):
    """Test that the module can write a simple stream input file"""

    # write test plt file to working dir
    confirm_name = 'test_origen_input.inp'
    confirm = open(general.datapath+confirm_name).read()
    test = origen_inst.make_file(general.wdpath, general.name, simple_stream,
                                 0.0, 5.0, 5.0, 'I', cooling=True)
    test = open(general.wdpath+'/'+test).read()
    tag = test.find('=origen') # need to omit header
    assert test[tag:] == confirm[tag:]


def test_run_origen(general, simple_stream, origen_inst):
    """Test that origen class can run scale"""

    test = origen_inst.make_file(general.wdpath, general.name, simple_stream,
                                 0.0, 5.0, 5.0, 'I', cooling=True)
    origen_inst.run_origen(general.wdpath, test)
    #assert '*.plt' in os.listdir(general.wdpath)
    assert len([x for x in os.listdir(general.wdpath)
                if test[:-4] in x and '.plt' in x]) == 2
    scale_msg = open(os.path.join(general.wdpath,test[:-4]+'.msg')).read()
    assert 'error' or 'ERROR' or 'Error' in scale_msg


def test_scale_path():
    """Test that scale has been added to path"""
    
    # will fail with OSError if not added to bash PATH
    path = os.environ['PATH'].split(os.pathsep)
    assert any('SCALE' in x for x in path), 'No path specified to SCALE'
    path_to_scale = [x for x in path if 'SCALE' in x]
    if(len(path_to_scale)>1):
        path_to_scale = [x for x in path_to_scale if 'bin' in x]
    assert len(path_to_scale) > 0, 'need better path management to SCALE dir'
    path_to_scale = os.path.join(path_to_scale[0], 'scalerte')
    try:
        assert os.path.exists(path_to_scale) # osx
    except:
        assert os.path.exists(path_to_scale+'.exe') # windows


def test_format_materials():
    """Test that the materials can be properly formatted"""

    mat = {'a':3.14, 'b':69.0, 'c':420.2, 'd':8008.5}
    x = origen.format_materials(mat)
    entries = x.split(' ')[1:len(x.split(' '))-1]
    assert len(entries) == len(mat)
    for i in range(0, len(entries)):
        assert float(entries[i][2:])== mat[entries[i][0]]


def test_format_time():
    """Test that the time input is properly formatted"""

    assert origen.format_time(100.0, 9, 'I') == '[8I 1 100.0]'
    assert origen.format_time(100.0, 9, 'L') == '[8L 1 100.0]'
    #assert origen.format_time(0.0, 0, 'L') == '0.0'


def test_format_opus_nuclides():
    """Test that the nuclides are properly formatted"""

    x = origen.format_opus_nuclides(origen.nuclides)
    assert len(origen.nuclides.split(' ')) == x[1]
    assert all([len(y) < 72 for y in x[0].split('\n')])


def test_origen_attributes():
    """Test the required attributes of the origen module"""
    
    assert hasattr(origen, 'preamble_temp')
    assert hasattr(origen, 'decay_case_temp')
    assert hasattr(origen, 'opus_temp')
    assert hasattr(origen, 'nuclides')


###############################################################################
