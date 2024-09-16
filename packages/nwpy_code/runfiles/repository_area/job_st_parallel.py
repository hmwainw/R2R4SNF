################################################################################
"""
    Last modified: September 27, 2018

    @author: Milos Atz <milos.atz@berkeley.edu>

    comments: Parallel run file for repository footprint surface storage
    time calculation with ipyparallel on Savio

"""
################################################################################
import numpy as np
import os
import imp
import types
import copy_reg
import ipyparallel as ipp
from footprint import repository
from footprint.repository import calculate_surface_storage_times
################################################################################


"""
Three strategies are available to ensure that repository thermal limits
are not broken by waste emplacement:

1) decrease package loading
2) increase surface storage time
3) increase spacing between packages

In the parametric analysis controlled by this run file, the surface storage
time required to (a) make disposal possible, and (b) reach the minimum 
repository footprint, is calculated.

"""


#-------------------------------------------------------------------------------
# WASTE TYPE
# The waste type controls the decay heat rate of the packages, which in turn
# affects the package spacing, possible loading options, and the required
# surface storage time.
#-------------------------------------------------------------------------------
waste_type = ['uox', 'mox', 'coex', 'nuex', 'ecc', 'ecm']


#-------------------------------------------------------------------------------
# REPOSITORY TYPE
# The repository host rock and disposal concept informs the heat transfer
# properties of the rock and engineered barrier system materials, as well as
# the thermal constraint that must be met.
#-------------------------------------------------------------------------------
rock_type = ['granite', 'clay', 'salt']


#-------------------------------------------------------------------------------
# PACKAGE LOADING
# Increasing the amount of waste loaded per package can decrease costs, but
# it also makes it more difficult to meet the thermal constraints. Only SNF
# is considered for variable package loading, as the pour-canister loading
# for HLW is already optimized for heat content (i.e. the decay heat of a
# glass waste form must be less than some value to prevent devitrification).
#-------------------------------------------------------------------------------
n_wf = {}
n_wf['uox'] = [1, 2, 3, 4, 12]
n_wf['mox'] = [1, 2, 3, 4, 12]
n_wf['coex'] = [1]
n_wf['nuex'] = [1]
n_wf['ecc'] = [1]
n_wf['ecm'] = [1]


#-------------------------------------------------------------------------------
# STORAGE TIME
# Waste can be allowed to cool in interim storage on the surface before
# emplacement in a repository. During this time, shorter-lived fission
# products can decay, making it easier to meet repository thermal constraints.
#-------------------------------------------------------------------------------
storage_time = [5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0,
                100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0,
                350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 700.0, 800.0, 900.0,
                1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0] # y


################################################################################


#-------------------------------------------------------------------------------
# The "pickle" module does not immediately pickle instance methods. In order
# to force it to do that, we need to change its behavior using copy_reg
#-------------------------------------------------------------------------------
def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)
#-------------------------------------------------------------------------------


client = ipp.Client()
lv = client.load_balanced_view()
tasks = []
for r in rock_type:
    c = repository.Array(r) # instantiate repository array
    for w in waste_type: # specify waste type...
        for n in n_wf[w]: # ... and package loading ...
            file = w+'-'+str(n)+'.csv' # load waste from file
            c.load_waste_from_file('./footprint/data/waste/'+file)
            tasks.append(lv.apply(calculate_surface_storage_times(c,r,w,n))


result = [task.get() for task in tasks]  # blocks until all results are back
print('done')


################################################################################


