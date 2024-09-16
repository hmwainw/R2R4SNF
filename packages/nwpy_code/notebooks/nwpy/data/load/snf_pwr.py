###############################################################################
"""
    Last modified on May 15, 2019
        
    @author: matz
        
    comment: data file containing waste form loading information for loading
             used pwr fuel assemblies into canisters
    
"""
###############################################################################
# SPENT FUEL: Pressurized Waste Reactor (PWR)
###############################################################################
hm_per_asm = 460000.0 # g/asm (460 kg)
asm_per_canister = 4 # DEFAULT VALUE; user can/should override
canister = {} # values in m
canister[1] = {'Diameter': 0.41, 'Length': 5.0, 'Thickness':0.125}
canister[2] = {'Diameter': 0.82, 'Length': 5.0, 'Thickness':0.125}
canister[3] = {'Diameter': 0.82, 'Length': 5.0, 'Thickness':0.125}
canister[4] = {'Diameter': 0.82, 'Length': 5.0, 'Thickness':0.125}
canister[12] = {'Diameter': 1.29, 'Length': 5.0, 'Thickness':0.125}
canister[21] = {'Diameter': 1.60, 'Length': 5.0, 'Thickness':0.125}
# hm per asm: ORNL/TM-9591/V1&R1, Table 3.1
# canister dimensions from SNL2012 repository thermal report
###############################################################################
