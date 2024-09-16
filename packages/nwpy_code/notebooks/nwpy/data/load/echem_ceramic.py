################################################################################
"""
    Last modified on May 16, 2019
        
    @author: matz
        
    comment: data file containing waste form loading information for ceramic
             waste form from electrochemical (echem) reprocessing.
    
"""
################################################################################
# CERAMIC (sodalite from pyroprocessing)
canister = {}
canister['Diameter'] = 0.6096 # m
canister['Length'] = 4.572 # m
canister['Mass limit'] = 2.9e6 # grams
canister['Thickness'] = 0.06 # 1 cm thick canister and 5 cm overpack
################################################################################
sodalite = {}
# fraction of fp in salt (no ion exchange or distillation)
sodalite['salt fp fraction'] = 0.10
# fraction of fp in salt (if ion exchange and distillation)
sodalite['salt high fp fraction'] = 0.95
# fraction of salt in salt-loaded zeolite
sodalite['SLZ salt fraction'] = 0.10
# after mixing with glass, ceramic waste form has fraction of zeolite
sodalite['CWF zeolite fraction'] = 0.75
# canister mass
sodalite['Canister mass limit'] = canister['Mass limit']
################################################################################