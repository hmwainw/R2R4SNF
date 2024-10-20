###############################################################################
"""
    Last modified on January 5, 2022
    
    @author: matz
    
    comment: data file containing waste package loading information 
             for TRISO-bearing fuel pebbles into waste packages
             
"""
###############################################################################
# SPENT FUEL: TRISO particles (LIFE FFH)
# TRISO particles are embedded into a 2 cm diameter carbon pebble with 15%
# packing fraction.
from __future__ import absolute_import
import fnmatch
import os
import imp
import numpy as np
# IMPORT PWR SNF PACKAGE DATA
dir_path = os.path.dirname(os.path.realpath(__file__))
matches = []
for root, dirnames, filenames in os.walk(dir_path):
    for filename in fnmatch.filter(filenames, 'snf_pwr.py'):
        matches.append(os.path.join(root, filename))
matches = [m for m in matches if not any([x in m for x in ['/old', '/test']])]
assert len(matches)==1
snf_pwr = imp.load_source('snf_pwr', matches[0])
###############################################################################
# Fuel pebble, triso particle characteristics
d_pebble = 2.0 # cm
d_triso = 0.105 # cm
d_kernel = 600.0/1e4 # cm
dens_kernel = 10.0 # g/cc
vfrac_triso = 0.15
# Calculate HM/pebble (g)
v_pebble = (4/3.0)*np.pi*(d_pebble/2)**3
v_triso = (4/3.0)*np.pi*(d_triso/2)**3
triso_per_pebble = v_pebble*vfrac_triso/v_triso
hm_per_triso = dens_kernel*((4/3.0)*np.pi*(d_kernel/2)**3) # g
hm_per_asm = hm_per_triso*triso_per_pebble
# Calculate pebbles per package; create canister dict
# ASSUME SAME PACKAGE DIMENSIONS AS 4-asm, 12-asm, and 21-asm PWR
vf = 0.4 # void fraction
canister = {}
for npwr in [1,4,12,21]:
    d_inner = (snf_pwr.canister[npwr]['Diameter']-
               snf_pwr.canister[npwr]['Thickness'] * 2)
    r_inner = d_inner/2.0
    l_inner = (snf_pwr.canister[npwr]['Length']-
               2*snf_pwr.canister[npwr]['Thickness'])
    v_inner = np.pi*r_inner**2*l_inner*1e6 # m3 -> cm3
    print(v_inner)
    pebbles_per_pkg = int(np.floor(v_inner*(1-vf)/v_pebble))
    canister[pebbles_per_pkg] = snf_pwr.canister[npwr]
# DEFAULT VALUE: 4 asm PWR package (should be second smallest)
asm_per_canister = sorted(canister.keys())[1] #max(canister.keys())
# References
# [1] https://escholarship.org/uc/item/6v84b2c6
# [2] http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.1018.3700&rep=rep1&type=pdf
###############################################################################
