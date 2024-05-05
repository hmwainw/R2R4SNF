###############################################################################
"""
    Last modified on May 16, 2019
        
    @author: matz
        
    comment: data file containing information for the loading of used HTGR
             fuel blocks/pebbles into waste canisters.
             
"""
###############################################################################
from __future__ import absolute_import
import numpy as np
###############################################################################
# SPENT FUEL: High Temperature Gas Reactor (HTGR)
# Assume Fort-Saint Vrain prismatic block-type fuel
###############################################################################
# DIRECT DISPOSAL OF PRISMATIC BLOCKS
hm_per_asm = 7200.0 # g initial HM per prismatic block (7.2 kg)
asm_per_canister = 18 # DEFAULT VALUE; user can/should override
canister = {}
canister[6] = {'Diameter': 0.66+0.125,
               'Length': 4.7625+0.125,
               'Thickness':0.125} # m
canister[18] = {'Diameter': 0.9298+0.125,
                'Length': 4.7625+0.125,
                'Thickness':0.125} # m
canister[42] = {'Diameter': 1.1891+0.125,
                'Length': 4.7625+0.125,
                'Thickness':0.125} # m
###############################################################################
# CONSOLIDATION OF FUEL ELEMENTS, REMOVAL OF GRAPHITE
# The block fuel channels are filled with fuel elements containing TRISO
# fuel particles embedded in graphite. Removal of the surrounding graphite
# to dispose of only the fuel elements increases package loading. None of
# the inserts that divide the hexagonal blocks are required, so more volume
# is free to be occupied by waste. Conservatively, the void volume caused
# by imperfect packing in the canister is assumed to be 40%. The "fuel_dens"
# (units: gHM/m3) multiplied by the canister inner volume yields the package
# loading.
v_ch = 9.846e-5 # m^3 per cylindrical fuel channel
n_ch = 210 # per block
vf = 0.4 # Void fraction loading elements into pkgs: 1-vf = v_fuel/v_can
fuel_dens = hm_per_asm*(1-vf)/(v_ch*n_ch)
###############################################################################
# Sources:
# HM Loading: https://inldigitallibrary.inl.gov/sites/sti/sti/4654912.pdf
# Canister and prismatic block dimensions: ORNL/TM-12027
###############################################################################