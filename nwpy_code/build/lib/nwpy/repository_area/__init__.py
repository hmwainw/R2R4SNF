###############################################################################
"""
    Last edited on June 28, 2019

    @author: matz

    comments: Python set up directory paths
    from https://stackoverflow.com/questions/4519127/
    setuptools-package-data-folder-location

    
"""
###############################################################################
import os
import nwpy.repository_area.repository as repository
###############################################################################
_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path, test=False):
    if(test==True): # direct the import to the testdata directory
        return(os.path.join(_ROOT, 'tests/testdata', path))
    else:
        return(os.path.join(_ROOT, 'data', path))
###############################################################################
