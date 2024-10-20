################################################################################
"""
    Last modified on May 20, 2019
        
    @author: matz
        
    comment: data file containing waste form loading information for loading
             melt refining crucible skull oxide slag into canisters

"""
################################################################################
# OXIDE WASTE (crucible skull) from melt-refining
# If melt refinining is to be used in a continuous recycle process, the skull
# would need to be processed to further recover U and Pu. However, because melt
# refining is only used in the limited-recycle options and its use is limited
# to 3 recycle steps, it is assumed that no skull processing is required.
# Therefore, I assume the skull is removed from the crucible via oxidation.
# Then, the resulting oxides are disposed of in a glass similar to UREX glass.
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
