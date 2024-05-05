###################################################################################
"""
    Last edited on July 11, 2018

    @author: matz

    comments: Test the methods and attributes of the FuelCycle class
    
"""
###################################################################################
import pytest
import os
import fuelcycle
from fuelcycle import fuelcycle
###################################################################################


@pytest.fixture
def eg01():
    print('\n'+'Setup EG01 fc instance')
    return(fuelcycle.FuelCycle('eg01', test=True))


###################################################################################


def test_init_invalid_evaluation_group():
    """Test that an nonexistant evaluation group results in proper error"""
    
    # actually testing _get_fc_data
    with pytest.raises(IOError):
        fuelcycle.FuelCycle('butt', 1, test=True) # no eg named 'butt'
        fuelcycle.FuelCycle('eg41', 1, test=True) # no eg 41


def test_init_paths(eg01):
    """Test that stage instance file paths all exist"""
    
    #assert os.path.exists(eg01.wdpath), 'path to working dir failed'
    assert os.path.exists(eg01.datapath),'path to data dir failed'


###################################################################################