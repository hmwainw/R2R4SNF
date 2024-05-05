################################################################################
"""
    Last modified on May 16, 2019
        
    @author: matz
        
    comment: data file containing waste form information for loading
             melt refining gas wastes (filters, etc) into canisters

"""
################################################################################
# GAS WASTE (Cs filters) from melt-refining
# Alkali, halide, and noble gas FP are volatilized from the fuel melt in the
# melt-refining process. The noble gases are held in tanks for decay (Kr) and
# controlled release. The alkali gases are trapped in filters as they leave
# the melt crucible. The halogens (notably iodine) are also captured. The gases
# with the highest activity by far are the alkali gases, namely Cs, and so
# those are considered here.
################################################################################
from __future__ import absolute_import
import numpy as np
################################################################################
canister = {}
canister['Diameter'] = 0.6096 # m
canister['Length'] = 4.572 # m
canister['Mass limit'] = 2.9e6 # grams
canister['Thickness'] = 0.06 # 1 cm thick canister, 5cm overpack
################################################################################
# cs loaded onto molecular sieve bed, beds crushed and put into canisters
# low loading value wrt temperature is assumed due to high process temperature
# for canister loading, volume is limiting
bed_loading = 3.14 # g alkali per cubic inch; molecular sieve
v_can = np.pi*(canister['Diameter']/2.0)**2*canister['Length'] # m^3
alkali_per_canister = bed_loading*1e6*v_can/16.3871
################################################################################
# 1 inch = 16.3871 cc
# https://www.tandfonline.com/doi/pdf/10.13182/NSE61-A25868?needAccess=true
################################################################################
