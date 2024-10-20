###############################################################################
"""
    Last modified on May 20, 2019
        
    @author: matz
        
    comment: data file containing waste form loading information for loading
             used sfr metal fuel assemblies into canisters
    
"""
###############################################################################
# SPENT FUEL: Sodium Fast Reactor (SFR metal fuel)
# Assume waste package is the same for ADS
###############################################################################
from __future__ import absolute_import
import fnmatch
import os
import imp
matches = []
for root, dirnames, filenames in os.walk('.'):
    for filename in fnmatch.filter(filenames, 'snf_ads.py'):
        matches.append(os.path.join(root, filename))
matches = [m for m in matches if not any([x in m for x in ['old', 'test']])]
assert len(matches)==1
canister = imp.load_source('snf_ads', matches[0]).canister
###############################################################################
hm_per_asm = 97700.0 # g/asm (from SNL 2011)
asm_per_canister = 7 # DEFAULT VALUE; user can, should override
# hm_per_asm = 13.1*1e6/180 # g/asm ~ 72.8 kg/asm (ABR-1000 reference)
# ABR-1000: https://publications.anl.gov/anlpubs/2017/04/134264.pdf
# hexagonal assembly: L ~477.52 cm; d_short ~15.71 cm; d_long ~18.14 cm
###############################################################################
