###############################################################################
"""
    Last edited on May 22, 2019

    @author: matz

    comments: Test the methods in the loading module
    
"""
###############################################################################
import pytest
import numpy as np
import pandas as pd
from copy import deepcopy
from nwpy import loading
from nwpy import stream
from nwpy import nuclidedata
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
    print('\n'+'Setup mock SNF Stream instance')
    c = pd.DataFrame({5.0: [1e5, 96e5, 1e5, 1e5, 1e6]},
                     index=['u235','u238','pu239','sr90','cs137'])
    h = pd.DataFrame({5.0: [1e3, 0.0, 1e3, 10e4, 10e4]},
                     index=['u235','u238','pu239','sr90','cs137'])
    return(stream.Stream(sum(c[5.0]), comp=c, heat=h, id='snf', form='mocksnf',
                         evaluationgroup='test', stagenumber=1))


@pytest.fixture
def simple_uranium(): # fake composition, fake heat
    print('\n'+'Setup simple mock HLW Stream instance')
    c = pd.DataFrame({5.0: [0.1, 1.0]}, index=['u235','u238'])
    h = pd.DataFrame({5.0: [0.0, 0.0]}, index=['u235','u238'])
    return(stream.Stream(1.1, comp=c, heat=h, id='hlw', form='glass',
                         evaluationgroup='test', stagenumber=1))


@pytest.fixture
def mock_hlw(): # fake composition, fake heat
    print('\n'+'Setup mock HLW Stream instance')
    c = pd.DataFrame({5.0: [0.01, 1.0, 0.01, 1.0, 1.0]},
                     index=['u235','u238','pu239','sr90','cs137'])
    h = pd.DataFrame({5.0: [0.0, 0.0, 0.01, 10.0, 10.0]},
                     index=['u235','u238','pu239','sr90','cs137'])
    return(stream.Stream(3.02, comp=c, id='hlw', heat=h, form='mockhlw',
                         evaluationgroup='test', stagenumber=1))


@pytest.fixture
def mock_metal_fp(): # fake composition, fake heat
    print('\n'+'Setup mock Stream instance of metal FP')
    c = pd.DataFrame({5.0: [2.0, 1.5, 0.3]}, index=['zr95', 'mo95', 'ru103'])
    h = pd.DataFrame({5.0: [11.0, 0.0, 3.2]}, index=['zr95', 'mo95', 'ru103'])
    return(stream.Stream(3.8, comp=c, heat=h, id='hlw', form='metal',
                         evaluationgroup='test', stagenumber=1))


@pytest.fixture
def mock_gas_fp(): # fake composition, fake heat
    print('\n'+'Setup mock Stream instance of gaseous FP')
    c = pd.DataFrame({5.0: [2.0, 1.5, 0.3]}, index=['cs137', 'cs135', 'rb85'])
    h = pd.DataFrame({5.0: [9.1, 0.05, 4.2]}, index=['cs137', 'cs135', 'rb85'])
    return(stream.Stream(3.8, comp=c, heat=h, id='hlw', form='gas',
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
    assert round(test.mass*test.number,2) == round(mock_snf.mass, 2)
    assert (round(sum(test.heat[5.0])*test.number,2) ==
            round(sum(mock_snf.heat[5.0]),2))


def test_general_variable_loading(mock_snf, loading_general):
    """Assert snf mass,heat equals canister mass,heat x canister number"""
    
    test1 = loading_general.load_waste(mock_snf, loading=1)
    assert test1.canister['Diameter'] == 0.41 # pwr single assembly canister
    assert round(test1.mass*test1.number,2) == round(mock_snf.mass, 2)
    assert (round(sum(test1.heat[5.0])*test1.number,2) ==
            round(sum(mock_snf.heat[5.0]),2))
    test2 = loading_general.load_waste(mock_snf, loading=2)
    assert round(test2.mass) == round(test1.mass*2)
    assert round(float(test2.number)) == round(test1.number/2.0)


def test_general_custom_can_and_loading(mock_snf, loading_general):
    """Assert snf mass,heat equals canister mass,heat x canister number"""

    c = {'Diameter':2.0, 'Length':4.5}
    asm = 5 # assemblies per canister
    tol = 0.1 # arbitrary
    test = loading_general.load_waste(mock_snf, loading=asm, can=c)
    confirm = 460000.0*asm # should be less than this
    assert confirm*(1-tol) <= test.mass
    assert test.canister == c


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


def test_oxidize(simple_uranium, loading_glass):
    """Ensure oxidize reports proper result"""
    
    hlw = deepcopy(simple_uranium)
    loaddata = loading_glass._get_loading_data(hlw.id, hlw.form)
    res_ox_dict = loading_glass._oxidize(hlw.comp[5.0], loaddata.oxide)
    assert len(res_ox_dict) == 1 # only one entry for uranium
    assert 'u3o8' in res_ox_dict.keys()
    m_uox = (hlw.comp[5.0]['u238']*842/3/238+
             hlw.comp[5.0]['u235']*833/3/235)
    assert round(m_uox,2)== round(res_ox_dict['u3o8'],2)
    

def test_glass_mass_and_heat(mock_hlw, loading_glass):
    """Ensure glass mass, heat equals canister mass,heat x canister number"""

    glass_hlw = deepcopy(mock_hlw)
    glass_hlw.form = 'glass'
    test = loading_glass.load_waste(glass_hlw)
    print(test.number)
    print(test.mass)
    print(sum(test.comp[5.0]))
    
    assert round(test.mass*test.number,2) == round(glass_hlw.mass,2)
    assert (round(sum(test.heat[5.0])*test.number,2) ==
            round(sum(glass_hlw.heat[5.0]),2))


def test_glass_user_loading(mock_hlw, loading_glass):
    """Ensure user can properly specify custom loading fraction"""
    
    glass_hlw = deepcopy(mock_hlw)
    glass_hlw.form = 'glass'
    t = loading_glass.load_waste(glass_hlw, loading_fraction=0.5)
    est = mock_hlw.mass/0.5/t.canister['Mass limit']
    assert np.ceil(est)==t.number
    # try messing it up
    with pytest.raises(AssertionError):
        test = loading_glass.load_waste(glass_hlw, loading_fraction=-0.5)
    with pytest.raises(AssertionError):
        test = loading_glass.load_waste(glass_hlw, loading_fraction=1.5)


###############################################################################
# METAL FROM ELECTROCHEMICAL REPROCESSING
###############################################################################


def test_ecmetal_instantiation(loading_echem_metal):
    """Make sure ec metal Loading instance has required methods on top of
    those from the general loading class"""

    assert hasattr(loading_echem_metal, '_import_constraints')
    assert hasattr(loading_echem_metal, '_update_constraints')
    assert hasattr(loading_echem_metal, '_evaluate_constraints')


def test_metal_mass_and_heat(mock_metal_fp, loading_echem_metal):
    """Ensure metal mass, heat equals canister mass,heat x canister number"""
    
    test = loading_echem_metal.load_waste(mock_metal_fp)
    assert round(test.mass*test.number,2) == round(mock_metal_fp.mass,2)
    assert (round(sum(test.heat[5.0])*test.number,2) ==
            round(sum(mock_metal_fp.heat[5.0]),2))


def test_metal_user_loading(mock_metal_fp, loading_echem_metal):
    """Ensure user can properly specify custom loading fraction"""
    
    t = loading_echem_metal.load_waste(mock_metal_fp, loading_fraction=0.5)
    est = mock_metal_fp.mass/0.5/t.canister['Mass limit']
    assert np.ceil(est)==t.number
    # try messing it up
    with pytest.raises(AssertionError):
        t=loading_echem_metal.load_waste(mock_metal_fp,loading_fraction=-0.5)
    with pytest.raises(AssertionError):
        t=loading_echem_metal.load_waste(mock_metal_fp,loading_fraction=1.5)


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
    test = loading_echem_ceramic.load_waste(ceramic_hlw)
    assert round(test.mass*test.number,2) == round(ceramic_hlw.mass,2)
    assert (round(sum(test.heat[5.0])*test.number,2) ==
            round(sum(ceramic_hlw.heat[5.0]),2))


def test_ceramic_user_loading(mock_hlw, loading_echem_ceramic):
    """Ensure user can properly specify custom loading fraction"""
    
    ceramic_hlw = deepcopy(mock_hlw)
    ceramic_hlw.form = 'ceramic'
    t = loading_echem_ceramic.load_waste(ceramic_hlw, loading_fraction=0.5)
    est = ceramic_hlw.mass/0.5/t.canister['Mass limit']
    assert np.ceil(est)==t.number
    # try messing it up
    with pytest.raises(AssertionError):
        t=loading_echem_ceramic.load_waste(ceramic_hlw,loading_fraction=-0.5)
    with pytest.raises(AssertionError):
        t=loading_echem_ceramic.load_waste(ceramic_hlw,loading_fraction=1.5)


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


def test_gas_mass_and_heat(mock_gas_fp, loading_meltrefining_gas):
    """Ensure ceramic mass,heat equals canister mass,heat x canister number"""

    mock_gas_fp.mass = sum(mock_gas_fp.comp[5.0])
    test = loading_meltrefining_gas.load_waste(mock_gas_fp, verbose=True)
    assert round(test.mass*np.ceil(test.number),2)==round(mock_gas_fp.mass,2)
    assert (round(sum(test.heat[5.0])*np.ceil(test.number),2) ==
            round(sum(mock_gas_fp.heat[5.0]),2))


def test_gas_user_loading(mock_gas_fp, loading_meltrefining_gas):
    """Ensure user can properly specify custom loading fraction"""
    
    t=loading_meltrefining_gas.load_waste(mock_gas_fp, loading_fraction=0.5)
    est = mock_gas_fp.mass/0.5/t.canister['Mass limit']
    assert np.ceil(est)==t.number
    # try messing it up
    with pytest.raises(AssertionError):
        t=loading_meltrefining_gas.load_waste(mock_gas_fp,loading_fraction=-0.5)
    with pytest.raises(AssertionError):
        t=loading_meltrefining_gas.load_waste(mock_gas_fp,loading_fraction=1.5)

###############################################################################
# SKULL FROM MELT REFINING
# Basically a copy of glass - need to make sure it loads the same methods
###############################################################################


def test_skull_instantiation(loading_meltrefining_skull, loading_glass):
    """Skull Loading instance should have same methods as AqGlass"""
    
    # ceramic loading only has a custom load_waste method
    assert all([x in dir(loading_glass) for x in dir(loading_meltrefining_skull)])


###############################################################################

