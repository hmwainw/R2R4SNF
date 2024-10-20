################################################################################
"""
    Last modified on May 16, 2019
        
    @author: matz
        
    comment: data file containing waste form loading information for 
             fluorapatite ceramic waste from reprocessing msr salt
    
"""
################################################################################
# CERAMIC (fluorapatite loaded from MSR salt)
canister = {}
canister['Diameter'] = 0.6096 # m
canister['Length'] = 4.572 # m
canister['Mass limit'] = 2.9e6 # grams
canister['Thickness'] = 0.06 # 1cm thick canister, 5cm overpack
################################################################################
# MSR SALT COMPOSITION (compositions are in mol fraction element fluorides)
salt = {}
salt['li'] = 0.717  # LiF
salt['be'] = 0.160  # BeF2
salt['th'] = 0.1205 # ThF4
salt['u'] = 0.0025  # UF4
molar_mass = 64.0485 # g/mol salt
salt_to_th = molar_mass*100/232/salt['th'] # g salt per gram th
################################################################################
fluorapatite = {}
# Max loading of FP in ceramic
fluorapatite['FP loading'] = 0.10
fluorapatite['Canister mass limit'] = canister['Mass limit']
################################################################################
