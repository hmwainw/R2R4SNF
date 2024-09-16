###############################################################################
"""
    Last edited on August 29, 2018

    @author: matz

    comments: Test the methods in the loading module
    
"""
###############################################################################
import pytest
import numpy as np
from copy import deepcopy
from nwpy.fuelcycle import loading
from nwpy.fuelcycle import stream
from nwpy.fuelcycle import nuclidedata
###############################################################################


@pytest.fixture
def loading_general(general):
    print('\n'+'Setup general loading instance')
    return(loading.Loading(data=general.data, datapath=general.datapath))


@pytest.fixture
def loading_glass(general_urex):
    print('\n'+'Setup UREX glass loading instance')
    return(loading.AqGlass(data=general_urex.data,
                           datapath=general_urex.datapath))


@pytest.fixture
def loading_echem_metal(dfbb):
    print('\n'+'Setup e-chem metal loading instance')
    return(loading.EcMetal(data=dfbb.data, datapath=dfbb.datapath))


@pytest.fixture
def loading_echem_ceramic(dfbb):
    print('\n'+'Setup e-chem ceramic loading instance')
    return(loading.EcCeramic(data=dfbb.data, datapath=dfbb.datapath))


@pytest.fixture
def loading_msr_ceramic(msr):
    print('\n'+'Setup msr ceramic loading instance')
    return(loading.MSRCeramic(data=msr.data, datapath=msr.datapath))


@pytest.fixture
def loading_msr_metal(msr):
    print('\n'+'Setup msr ceramic loading instance')
    return(loading.MSRMetal(data=msr.data, datapath=msr.datapath))


@pytest.fixture
def loading_meltrefining_gas(bnbsfr):
    print('\n'+'Setup melt refining gas loading instance')
    return(loading.CapturedCs(data=bnbsfr.data, datapath=bnbsfr.datapath))


@pytest.fixture
def loading_meltrefining_skull(bnbsfr):
    print('\n'+'Setup melt refining skull loading instance')
    return(loading.Skull(data=bnbsfr.data, datapath=bnbsfr.datapath))


@pytest.fixture
def mock_snf(): # fake composition, fake heat
    print('\n'+'Setup simple Stream instance')
    c = {'u235':1.0, 'u238':96.0, 'pu239':1.0,
        'sr90':1.0, 'cs137':1.0}
    h = {'u235':0.1, 'u238':0.0, 'pu239':1.0,
        'sr90':10.0, 'cs137':10.0}
    return(stream.Stream(100.0, c, id='snf', heat=h, form='mocksnf', time=5.0,
                         evaluationgroup='test', stagenumber=1))


@pytest.fixture
def mock_hlw(): # fake composition, fake heat
    print('\n'+'Setup simple Stream instance')
    c = {'u235':0.01, 'u238':1.0, 'pu239':0.01,
        'sr90':1.0, 'cs137':1.0}
    h = {'u235':0.0, 'u238':0.0, 'pu239':0.01,
        'sr90':10.0, 'cs137':10.0}
    return(stream.Stream(3.2, c, id='hlw', heat=h, form='mockhlw', time=5.0,
                         evaluationgroup='test', stagenumber=1))


###############################################################################
# USED NUCLEAR FUEL
###############################################################################


def test_general_instantiation(loading_general):
    """Make sure Loading instance has required methods"""
    
    assert hasattr(loading_general, 'load_waste')
    assert hasattr(loading_general, '_get_loading_data')
    assert hasattr(loading_general, '_make_wasteform')


def test_general_mass_and_heat(mock_snf, loading_general):
    """Assert snf mass,heat equals canister mass,heat x canister number"""

    test = loading_general.load_waste(mock_snf)
    assert round(test.mass*test.number,2) == round(mock_snf.mass,2)
    assert (round(sum(test.heat.values())*test.number,2) ==
            round(sum(mock_snf.heat.values()),2))


###############################################################################
# GLASS FROM AQUEOUS REPROCESSING
###############################################################################


def test_aqglass_instantiation(loading_glass):
    """Make sure AqGlass Loading instance has required methods on top of 
    those from the general loading class"""

    assert hasattr(loading_glass, '_import_constraints')
    assert hasattr(loading_glass, '_update_constraints')
    assert hasattr(loading_glass, '_evaluate_constraints')
    assert hasattr(loading_glass, '_oxidize')


def test_oxidize(mock_hlw, loading_glass):
    """Ensure oxidize reports proper result"""
    
    glass_hlw = deepcopy(mock_hlw)
    glass_hlw.form = 'glass'
    glass_hlw.comp = {'u235':0.1, 'u238':1.0} # real simple
    loaddata = loading_glass._get_loading_data(glass_hlw.id, glass_hlw.form)
    res_ox_dict = loading_glass._oxidize(glass_hlw.comp, loaddata.oxide)
    assert len(res_ox_dict) == 1
    assert 'u3o8' in res_ox_dict.keys()
    m_uox = glass_hlw.comp['u238']*842/3/238+glass_hlw.comp['u235']*833/3/235
    assert round(m_uox,2)== round(res_ox_dict['u3o8'],2)
    

def test_glass_mass_and_heat(mock_hlw, loading_glass):
    """Ensure glass mass, heat equals canister mass,heat x canister number"""

    glass_hlw = deepcopy(mock_hlw)
    glass_hlw.form = 'glass'
    test = loading_glass.load_waste(glass_hlw)
    assert round(test.mass*test.number,2) == round(glass_hlw.mass,2)
    assert (round(sum(test.heat.values())*test.number,2) ==
            round(sum(glass_hlw.heat.values()),2))


###############################################################################
# METAL FROM ELECTROCHEMICAL REPROCESSING
###############################################################################


def test_ecmetal_instantiation(loading_echem_metal):
    """Make sure ec metal Loading instance has required methods on top of
    those from the general loading class"""

    assert hasattr(loading_echem_metal, '_import_constraints')
    assert hasattr(loading_echem_metal, '_update_constraints')
    assert hasattr(loading_echem_metal, '_evaluate_constraints')


def test_metal_mass_and_heat(mock_hlw, loading_echem_metal):
    """Ensure metal mass, heat equals canister mass,heat x canister number"""
    
    metal_hlw = deepcopy(mock_hlw)
    metal_hlw.form = 'metal'
    metal_hlw.comp = {'zr95': 2.0, 'mo95':1.5, 'ru103':0.3}
    metal_hlw.heat = {'zr95': 11.0, 'mo95':0.0, 'ru103':3.2}
    metal_hlw.mass = sum(metal_hlw.comp.values())
    test = loading_echem_metal.load_waste(metal_hlw)
    assert round(test.mass*test.number,2) == round(metal_hlw.mass,2)
    assert (round(sum(test.heat.values())*test.number,2) ==
            round(sum(metal_hlw.heat.values()),2))


###############################################################################
# CERAMIC FROM ELECTROCHEMICAL REPROCESSING
###############################################################################


def test_ecceramic_instantiation(loading_echem_ceramic):
    """Make sure ec ceramic Loading instance has required methods on top of
        those from the general loading class"""
    
    # ceramic loading only has a custom load_waste method
    assert hasattr(loading_echem_ceramic, 'load_waste')


def test_ceramic_mass_and_heat(mock_hlw, loading_echem_ceramic):
    """Ensure ceramic mass,heat equals canister mass,heat x canister number"""
    
    ceramic_hlw = deepcopy(mock_hlw)
    ceramic_hlw.form = 'ceramic'
    ceramic_hlw.comp = {'u238': 0.1, 'cs135':1.5, 'gd154':0.3}
    ceramic_hlw.heat = {'u238': 0.0, 'cs135':9.1, 'gd154':1.2}
    ceramic_hlw.mass = sum(ceramic_hlw.comp.values())
    test = loading_echem_ceramic.load_waste(ceramic_hlw)
    assert round(test.mass*test.number,2) == round(ceramic_hlw.mass,2)
    assert (round(sum(test.heat.values())*test.number,2) ==
            round(sum(ceramic_hlw.heat.values()),2))


###############################################################################
# CERAMIC FROM MSR REPROCESSING
###############################################################################

###############################################################################
# WASTES FROM MSR SALT TREATMENT
###############################################################################


###############################################################################
# CS/RB OFF GAS FROM MELT REFINING
###############################################################################


def test_csrb_instantiation(loading_meltrefining_gas):
    """Make sure cs-rb Loading instance has required methods on top of
        those from the general loading class"""
    
    # melt refining gas loading only has a custom load_waste method
    assert hasattr(loading_meltrefining_gas, 'load_waste')


def test_gas_mass_and_heat(mock_hlw, loading_meltrefining_gas):
    """Ensure ceramic mass,heat equals canister mass,heat x canister number"""

    gas_hlw = deepcopy(mock_hlw)
    gas_hlw.form = 'gas'
    gas_hlw.comp = {'cs137': 2.0, 'cs135':1.5, 'rb85':0.3}
    gas_hlw.heat = {'cs137': 9.1, 'cs135':0.05, 'rb85':4.2}
    gas_hlw.mass = sum(gas_hlw.comp.values())
    test = loading_meltrefining_gas.load_waste(gas_hlw)
    assert round(test.mass*np.ceil(test.number),2) == round(gas_hlw.mass,2)
    assert (round(sum(test.heat.values())*np.ceil(test.number),2) ==
            round(sum(gas_hlw.heat.values()),2))


###############################################################################
# SKULL FROM MELT REFINING
# Basically a copy of glass - need to make sure it loads the same methods
###############################################################################


def test_skull_instantiation(loading_meltrefining_skull, loading_glass):
    """Skull Loading instance should have same methods as AqGlass"""
    
    # ceramic loading only has a custom load_waste method
    assert all([x in dir(loading_glass) for x in dir(loading_meltrefining_skull)])


###############################################################################

