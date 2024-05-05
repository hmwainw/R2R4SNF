###############################################################################
"""
    Last edited on February 19, 2019

    @author: matz

    comments: Test the methods in the nuclidedata module
    
"""
###############################################################################
import random
import string
import pytest
import numpy as np
from nwpy import nuclidedata
###############################################################################


def test_nuclide_info_id_length():
    """Test that determine_nuclide_info returns error if input is invalid"""
    
    with pytest.raises(AssertionError):
        nuclidedata.determine_nuclide_info('HE-4000') # Input too long


def test_nuclide_id_form():
    """Test that determine_nuclide_info parses random nuclide IDs"""
    
    for bupkis in range(0, 100):
        el = ''
        for i in range(int(np.random.random()*2+1)):
            el+=random.choice(string.ascii_letters.lower())
        if(np.random.random()>0.5):
            hyphen = '-'
        else:
            hyphen = ''
        mass = str(int(np.random.random()*100))
        if(np.random.random()>0.5):
            meta = 'M'
            mstatus = 'meta'
        else:
            meta = ''
            mstatus = 'not meta'
        id = el+hyphen+mass+meta
        assert nuclidedata.determine_nuclide_info(id)==(el,str(mass),mstatus)


def test_is_actinide():
    """Test for all actinides, show that no other nuclides are included"""
    
    for el in nuclidedata.Z.keys():
        if(nuclidedata.is_actinide(el)):
            assert(el in nuclidedata.actinides)
        else:
            assert(el not in nuclidedata.actinides)


def test_is_transuranic():
    """Test for all tru, show that no other nuclides are included"""
    
    for el in nuclidedata.Z.keys():
        if(nuclidedata.is_transuranic(el)):
            assert(el in nuclidedata.transuranics)
        else:
            assert(el not in nuclidedata.transuranics)


def test_is_minor_actinide():
    """Test for minor actinides, show that no other nuclides are included"""
    
    for el in nuclidedata.Z.keys():
        if(nuclidedata.is_minoractinide(el)):
            assert(el in nuclidedata.minor_actinides)
        else:
            assert(el not in nuclidedata.minor_actinides)


def test_nuclidedata_attributes():
    """Test that nuclidedata module has required attributes"""

    assert hasattr(nuclidedata, 'actinides')
    assert hasattr(nuclidedata, 'minor_actinides')
    assert hasattr(nuclidedata, 'transuranics')
    assert hasattr(nuclidedata, 'u3')
    assert hasattr(nuclidedata, 'Z')
    assert hasattr(nuclidedata, 'group_nuclides')


def test_isgroup():
    """Test that nuclide data function is_group properly sorts nuclides"""
    
    assert nuclidedata.is_group('co60', 'fp')
    assert not nuclidedata.is_group('PU-239', 'fp')
    assert nuclidedata.is_group('Ac', 'actinide')
    assert not nuclidedata.is_group('He-4', 'actinide')
    assert nuclidedata.is_group('AM241', 'tru')
    assert not nuclidedata.is_group('u-233', 'transuranic')
    assert nuclidedata.is_group('np237', 'ma')
    assert not nuclidedata.is_group('Ac', 'minor actinide')
    assert nuclidedata.is_group('u235', 'u3')
    assert not nuclidedata.is_group('th232', 'u3')

###############################################################################
