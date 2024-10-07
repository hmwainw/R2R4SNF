################################################################################
"""
    Last modified on May 16, 2019
        
    @author: matz
        
    comment: data file containing waste form loading information for metal
             waste from electrochemical (echem) reprocessing.
    
"""
################################################################################
from __future__ import absolute_import
import numpy as np
################################################################################
# METAL from electrochemical processing
canister = {}
canister['Diameter'] = 0.6096 # m
canister['Length'] = 3.048 # m
canister['Mass limit'] = 3.6e6 # grams
canister['Thickness'] = 0.06 # 1 cm thick canister, 5 cm overpack
################################################################################
# Assumptions:
# 1. Metal fuel has 10 wt% Zr; 0.1 = M_Zr/(M_Zr+M_HM)
zr_to_hm = 1/9.0
# 2. Cladding: assume all assemblies are the same
ht9_to_hm = 54.93/97.7 # kg clad / kgHm
# Inequality constraints
# 1. Nonzero waste mass         M_waste > 0
# 2. Nonzero fill (SS/Zr) mass  M_filler > 0
# 3. Mass limit: 3600 kg        M_waste + M_filler <= 3600 kg
# X. Minimum loading: 0.005     M_waste / (M_waste + M_filler) >= 0.001
metal = {}
metal['ineq'] = np.zeros((3,3))
metal['ineq'][0] = [-1.0,      0.0,   0.0]
metal['ineq'][1] = [0.0,       -1.0,  0.0]
metal['ineq'][2] = [1.0,       1.0,   canister['Mass limit']]
#metal['ineq'][3] = [0.005-1.0, 0.005, 0.0]
# Equality constraints
# 1. Waste form composition:
#   If HT-9 clad/hardware, waste form is SS-15Zr: Zr = 15 wt.%
#   NOTE: If Zry clad/hardware, wasteform is Zr-8SS: SS (use Fe+Cr+Ni) = 8 wt.%
metal['eq'] = np.zeros((1,3))
metal['eq'][0] = [0.0, 0.0, 0.0]
#
# KEY
metal['key'] = []
metal['key'].append('Nonzero waste mass')
metal['key'].append('Nonzero metal fill mass')
metal['key'].append('Max waste mass')
#metal['key'].append('Min waste loading')
metal['key'].append('Waste form composition')
################################################################################