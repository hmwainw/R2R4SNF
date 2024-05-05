###############################################################################
"""
    Last modified on May 15, 2019
        
    @author: matz
        
    comment: data file containing information for loading used hwr fuel 
             assemblies into waste canisters.
             
"""
###############################################################################
# SPENT FUEL: Heavy Water Reactor (HWR)
# Value                 Canada      S. Korea
# --------------------------------------------
# Outer Diamater (m)    1.168       1.12
# Length (m)            3.867       4.83
# Fuel bundles          324         297
# HM/bundle (kg)        NA          ~19.0
###############################################################################
hm_per_asm = 20000.0 # g/bundle (20 kg)
asm_per_canister = 324 # DEFAULT VALUE; user can/should override
canister = {} # dimensions in meters
canister[108] = {'Diameter': 1.168, 'Length': 1.7962, 'Thickness':0.1}
canister[216] = {'Diameter': 1.168, 'Length': 2.8316, 'Thickness':0.1}
#canister[297] = {'Diameter': 1.12, 'Length': 4.83, 'Thickness':0.125} # S. KOREA
canister[324] = {'Diameter': 1.168, 'Length': 3.867, 'Thickness':0.1} # CANADA
canister[432] = {'Diameter': 1.168, 'Length': 4.9024, 'Thickness':0.1}
canister[540] = {'Diameter': 1.168, 'Length': 5.9378, 'Thickness':0.1}
# Sources:
# Korea: https://www.tandfonline.com/doi/pdf/10.1080/18811248.2007.9711407
# Korea (2): https://www.nap.edu/read/11320/chapter/7#56 (see Fig 3)
# Canada: https://www.nwmo.ca/~/media/Site/Files/PDFs/2015/11/09/12/54/662_6-7StatusofStorageDisposalandTransportationContainersfortheManagementofUsedNuclearFuel.ashx?la=en
# Canada hypothetical canister dimensions are calculated assuming 1, 2, 4, or
# 5 108-bundle baskets are stacked in a canister. Each basket contains two
# layers of 54 bundles. The bundles are 495 mm in length and the basket is
# 1035.4 mm in length. In the 3-basket package, there is 760.8 mm of extra,
# non-fuel material that adds to the length; assuming all thicknesses stay the
# same, this value is added to the height of the stacked baskets to determine
# the canister height.
# basket height: https://www.nwmo.ca/~/media/Site/Files/PDFs/2015/11/17/23/25/1002_deepgeologicrepository_figures.ashx?la=en
###############################################################################
