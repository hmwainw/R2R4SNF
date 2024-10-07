###############################################################################
"""
    Last edited on February 19, 2019
    
    @author: matz

    comments: Test fixtures (stages) to be used to test the fuelcycle code
    
"""
###############################################################################
import pytest
from nwpy import stage
###############################################################################


@pytest.fixture
def general():
    """Reactor with general-form isotopics and no separation"""
    
    print('\n'+'Setup general reactor instance')
    name = 'rx-gen_sep-none'
    s = stage.Stage(name, 1, test=True)
    return(s)



@pytest.fixture
def general_urex():
    """General case, two stages, UREX reprocessing"""
    
    print('\n'+'Setup general reactor instance (UREX)')
    name = 'rx-gen_sep-urex'
    s = stage.Stage(name, 1, test=True)
    return(s)


@pytest.fixture
def dfbb():
    """Reactor with df-bb-form isotopics and electrochemical separation"""

    print('\n'+'Setup DF-BB reactor instance')
    name = 'rx-dfbb_sep-echem'
    s = stage.Stage(name, 1, test=True)
    return(s)


@pytest.fixture
def ffh():
    """Reactor with ffh-form isotopics and no separation"""

    print('\n'+'Setup FFH reactor instance')
    name = 'rx-ffh_sep-none'
    s = stage.Stage(name, 1, test=True)
    return(s)


@pytest.fixture
def msr():
    """Reactor with msr-form isotopics and msr separations"""

    print('\n'+'Setup MSR reactor instance')
    name = 'rx-msr_sep-msr'
    s = stage.Stage(name, 1, test=True)
    return(s)


@pytest.fixture
def bnbsfr():
    """Reactor with bnbsfr-form isotopics and melt-refining separation"""

    print('\n'+'Setup BnB SFR reactor instance')
    name = 'rx-bnbsfr_sep-mr'
    s = stage.Stage(name, 1, test=True)
    return(s)


###############################################################################

