################################################################################
"""
    Last modified on May 20, 2019
    
    @author: matz
    
    comment: data file containing waste form loading information for loading
    thorex glass HLW into canisters
    
    """
################################################################################
# THOREX GLASS HLW
# Same waste form as for UREX process.
################################################################################
from __future__ import absolute_import
import fnmatch
import os
import imp
matches = []
for root, dirnames, filenames in os.walk('.'):
    for filename in fnmatch.filter(filenames, 'urex_glass.py'):
        matches.append(os.path.join(root, filename))
matches = [m for m in matches if not any([x in m for x in ['old', 'test']])]
assert len(matches)==1
urex_glass = imp.load_source('urex_glass', matches[0])
canister = urex_glass.canister
glass = urex_glass.glass
oxide = urex_glass.oxide
################################################################################