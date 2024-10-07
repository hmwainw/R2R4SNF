###############################################################################
"""
    Last modified on May 20, 2019
        
    @author: matz
        
    comment: Data file containing waste package loading information for
             used Pb-Bi-cooled ADS fuel assemblies into canisters
    
"""
###############################################################################
# Accelerator-driven system cooled by lead-bismuth-eutectic (LBE)
# This is assumed to apply to EG07 and EG16
#
# LBE ATW:
# from [1]: HM (BOEC): 2192 kg; 192 assemblies
# from [2]: HM (BOEC): 2256 kg; 192 assemblies
# from [1]:
#217 pins per assembly
#pitch-to-diameter ratio 1.727
#pin diameter 0.635 cm
#cladding thickness 0.056 cm
# hexagonal assemblies; d_long ~ 18.6 cm; d_short ~ 15.6 cm
#active height 106.68 cm (same as oxide ABR-1000);
#assume fuel assembly height is the same as well (477 cm)
# Sodium-ATW
# (Sodium substitute): HM (BOEC): 2899 kg; 192 assemblies [2]
# (Compact sodium):    HM (BOEC): 3373 kg; 192 assemblies [2]
#                      HM (BOEC): 2602 kg; 132 assemblies [2]
#                      HM (BOEC): 2320 kg; 120 assemblies [2]
###############################################################################
hm_per_asm = 2224.0*1e3/192 # ~2220 kg in the core over 192 asm; ~11.5 kg/asm
asm_per_canister = 7 # DEFAULT VALUE; user can/should override
canister = {} # values in m
canister[1] = {'Diameter': 0.46,'Length': 5.0,'Thickness':0.125} # d_in=0.21
canister[3] = {'Diameter': 0.66,'Length': 5.0,'Thickness':0.125} # d_in=0.41
canister[7] = {'Diameter': 0.77,'Length': 5.0,'Thickness':0.125} # d_in=0.52
canister[19] = {'Diameter': 1.11,'Length': 5.0,'Thickness':0.125} # d_in=0.86
# [1] https://www.tandfonline.com/doi/pdf/10.13182/NT135-162?needAccess=true
# [2] https://publications.anl.gov/anlpubs/2001/08/40152.pdf
# Waste package dimensions modified based on ORNL/TM-12027 (HTGR)
###############################################################################
