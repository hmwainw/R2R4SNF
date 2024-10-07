###############################################################################
"""
    Last edited on May 6, 2019

    @author: matz

    comments: Test fixtures to test the footprint code
    
"""
###############################################################################
import pytest
import os
import pandas as pd
from nwpy import stream
from nwpy import stage
from nwpy.repository_area import repository
###############################################################################


@pytest.fixture
def lwr_snf():
    """LWR SNF from WasteForm instance"""
    
    print('\n'+'Setup LWR SNF waste stream')
    s = stage.Stage('rx-gen_sep-none', 1, test=True)
    s.wdpath = os.path.join(s.wdpath[:s.wdpath.find('test')+4], 'footprint',
                            s.wdpath[s.wdpath.find('test')+5:])
    if(not os.path.exists(s.wdpath)):
        os.makedirs(s.wdpath)
    return(s.discharge_all_wastes(endtime=1e3, steps=10))


@pytest.fixture
def const_heat_wf():
    """Fake WF with constant heat"""
    
    print('\n'+'Setup fake constant-heat waste form')
    t = ['0.0', '1.0', '5.0', '10.0', '50.0', '100.0', '250.0', '500.0']
    h = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
    wf = stream.WasteForm(10.0, number=10.0, loading=4,
                          canister={'Diameter':0.96, 'Length':5.0},
                          composition={}, form='fake snf',
                          evaluationgroup='test', stagenumber=1)
    wf.heat = pd.DataFrame(data=[h], columns=t)
    wf.id = 'snf'
    return(wf)


@pytest.fixture
def granite_const_heat(const_heat_wf):
    """Granite repository with wf that has constant heat"""
    
    print('\n'+'Setup granite repository instance (constant heat)')
    return(repository.Repository('granite', const_heat_wf, st=0.0))


@pytest.fixture
def granite_no_waste():
    """Granite repository with wf that has constant heat"""
    
    print('\n'+'Setup granite repository instance (constant heat)')
    return(repository.Repository('granite', st=0.0))


###############################################################################
