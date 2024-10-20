###############################################################################
"""
    Last edited on May 10, 2019

    @author: matz

    comments: Test waste loading functions
    
"""
###############################################################################
import os
from nwpy.repository_area import waste
from nwpy.repository_area.__init__ import _ROOT
###############################################################################


datapath = os.path.join(_ROOT, 'tests', 'testdata')


def test_datafile_import():
    """Test the import on an example of the CSV data 
    compiled from the SNL 2011 reports"""
    
    filename = 'uox-4.csv'
    w = waste.Waste(os.path.join(datapath, filename))
    for x in ['id', 'n_wf', 'pkg', 'decay_heat']: # required
        assert x in dir(w)
    # make sure the decay heat result works - expect decreasing numbers
    assert w.decay_heat(5.0) - w.decay_heat(6.0) > 0.0


def test_wf_import(lwr_snf):
    """Generate data using nwpy and import the wasteform"""
    
    w = waste.Waste(lwr_snf)
    assert all(x in dir(w) for x in ['id', 'n_wf', 'pkg', 'decay_heat',
                                     'evaluationgroup', 'stagenumber'])
    # make sure the decay heat result works - expect decreasing numbers
    assert w.decay_heat(5.0) - w.decay_heat(6.0) > 0.0


def test_fake_wf(const_heat_wf):
    """Ensure that fake waste form can be imported"""

    w = waste.Waste(const_heat_wf)
    assert all(x in dir(w) for x in ['id', 'n_wf', 'pkg', 'decay_heat'])
    # make sure the decay heat result works - expect decreasing numbers
    assert w.decay_heat(5.0)==w.decay_heat(6.0)


def test_wf_file_import():
    """Read CSV file written by nwpy"""
    
    filename = 'eg01_file_io_test.csv'
    w = waste.Waste(os.path.join(datapath, filename))
    print(dir(w))
    # required information
    for x in ['id', 'n_wf', 'pkg', 'decay_heat']:
        assert x in dir(w)
    # optional information that I know is in the datafile
    for x in ['evaluationgroup', 'stagenumber']:
        assert x in dir(w)
    # make sure the decay heat result works - expect decreasing numbers
    assert w.decay_heat(5.0) - w.decay_heat(6.0) > 0.0


###############################################################################
