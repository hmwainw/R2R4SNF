###############################################################################
"""
    Last edited on June 28, 2019

    @author: matz

    comments: Package setup
    
"""
###############################################################################
import os
from nwpy import nuclidedata
from nwpy import stage
from nwpy import fuelcycle
from nwpy import repository_area
###############################################################################
_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path, test=False):
    if(test==True): # direct the import to the testdata directory
        return(os.path.join(_ROOT, 'tests/testdata', path))
    else:
        return(os.path.join(_ROOT, 'data', path))
###############################################################################