################################################################################
"""
    Last modified on May 16, 2019
        
    @author: matz
        
    comment: data file containing waste form loading information for 
             metal waste holding fp from msr salt treatment
    
"""
################################################################################
# CERAMIC (fluorapatite loaded from MSR salt)
canister = {}
canister['Diameter'] = 0.6096 # m
canister['Length'] = 4.572 # m
canister['Mass limit'] = 2.9e6 # grams
canister['Thickness'] = 0.06 # 1cm thick canister, 5 cm overpack
################################################################################
metal = {}
# Max loading of FP in metal
metal['FP loading'] = 0.10
metal['Canister mass limit'] = canister['Mass limit']
################################################################################
