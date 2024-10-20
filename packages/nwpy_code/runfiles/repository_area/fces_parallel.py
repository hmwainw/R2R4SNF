###############################################################################
"""
    Last modified: May 22, 2019
    
    @author: Milos Atz <milos.atz@berkeley.edu>
    
    comments: Parallel run file for repository footprint calculation
    with ipyparallel on Savio
    
"""
###############################################################################
import numpy as np
import os
import imp
import types
import numbers
import copy_reg
import ipyparallel as ipp
from collections import Sequence
import nwpy
from nwpy import fuelcycle
from nwpy import stage
from nwpy import repository_area
from nwpy.repository_area.repository import fces_app
from nwpy.repository_area.repository import fces_min_st
###############################################################################


"""
Three strategies are available to ensure that repository thermal limits
are not broken by waste emplacement:

1) decrease package loading
2) increase surface storage time
3) increase spacing between packages

In the parametric analysis controlled by this run file, the spacing between
packages is calculated for different combinations of waste, repository
concept, package loading, and surface storage time.

"""


###############################################################################
# The "pickle" module does not immediately pickle instance methods. In order
# to force it to do that, we need to change its behavior using copy_reg
#------------------------------------------------------------------------------
def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)


#------------------------------------------------------------------------------
# FUEL CYCLES
#------------------------------------------------------------------------------
fuelcycles = ['eg01','eg02','eg03','eg04','eg05','eg06','eg07','eg08','eg09','eg10',
              'eg11','eg12','eg13','eg14','eg15','eg16','eg17','eg18','eg19','eg20',
              'eg21','eg22','eg23','eg24','eg25','eg26','eg27','eg28','eg29','eg30',
              'eg31','eg32','eg33','eg34','eg35','eg36','eg37','eg38','eg39','eg40']
# ONCE-THROUGH
once_through = ['eg01','eg02','eg03','eg04','eg05','eg06','eg07','eg08']


#------------------------------------------------------------------------------
# REPOSITORY TYPE
# The repository host rock and disposal concept informs the heat transfer
# properties of the rock and engineered barrier system materials, as well as
# the thermal constraint that must be met.
#------------------------------------------------------------------------------
rocktypes = ['granite', 'clay', 'salt']
#rocktypes = ['granite']


#------------------------------------------------------------------------------
# PACKAGE LOADING
# Increasing the amount of waste loaded per package can decrease costs, but
# it also makes it more difficult to meet the thermal constraints. Only SNF
# is considered for variable package loading, as the pour-canister loading
# for HLW is already optimized for heat content (i.e. the decay heat of a
# glass waste form must be less than some value to prevent devitrification).
#------------------------------------------------------------------------------
def get_pkg_loading_options(stg, strm):
    """Get possible options for package loading based on waste
    discharged from the instantiated stage."""
    
    if(strm.id == 'hlw'):
        return([1])
    else: # snf
        x = nwpy.loading.Loading(datapath=stg.datapath,data=stg.data)
        dat = x._get_loading_data(strm.id, strm.form).canister.keys()
        if(all([i in dat for i in ['Mass limit', 'Diameter', 'Length']])):
            assert not any([isinstance(i, numbers.Number) for i in dat])
            return([1])
        else:
            return(dat)


def treat_like_msr(stg):
    """Determine whether the streams should be reprocessed and loaded
    like MSR streams, which require special treatment in this loop"""
    
    if(stg.data['reactor'] in ['msr', 'ffh'] and
       'salttreatment' in stg.data.keys()):
        return(True)
       else:
           return(False)


#------------------------------------------------------------------------------
# STORAGE TIME
# Waste can be allowed to cool in interim storage on the surface before
# emplacement in a repository. During this time, shorter-lived fission
# products can decay, making it easier to meet repository thermal constraints.
#
# The data produced by nwpy starts at 5 years of surface storage, because it
# is assumed that all waste is cooled for 5 years before further processing.
#------------------------------------------------------------------------------
cooling_time = 5.0 # years
storage_time = [0.0, 5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0, 95.0,
                120.0, 145.0, 170.0, 195.0, 220.0, 245.0, 270.0, 295.0, 345.0,
                395.0, 445.0, 495.0] # years
#storage_time = [95.0]#[15.0, 45.0, 95.0] # years


###############################################################################


client = ipp.Client()
lv = client.load_balanced_view()
tasks = []
for eg in fuelcycles:
    fc = fuelcycle.FuelCycle(eg)
    n_stg = fc.totalstages
    for stg in xrange(1, n_stg+1):
        s = stage.Stage(eg, stg)
        strm = s.discharge_streams()
        strm = s.cool(strm)#, rerun=False)
        strm = s.reprocess(strm)
        # waste loading - MSR is special case
        if(treat_like_msr(s)):
            strm = s.load_waste(strm)
        if(not isinstance(strm, Sequence)):
            strm = [strm]
        for w in strm:
            for pl in get_pkg_loading_options(s, w):
                if(not treat_like_msr(s)): # need to load the waste
                    wf = s.load_waste(w, loading=pl)
                else: # w is already loaded waste form
                    wf = w
                wf = s.decay(wf, endtime=1e4, steps=25)#, rerun=False)
                if(not isinstance(wf, Sequence)):
                    wf = [wf]
                    for wfi in wf:
                        for r in rocktypes:
                            c = repository_area.repository.Repository(r, wfi)
                            for st in storage_time:
                                c.update_st(st)
                                # I/O DATA: EG, STG, WF, LOAD, N_PKGS, ROCK, ST
                                csvdata = [eg, str(stg), wfi.form, str(pl),
                                           str(wfi.number), r, str(st+cooling_time)]
                                print(csvdata)
                                tasks.append(lv.apply(fces_min_st,c,csvdata))
                                if(eg not in once_through):
                                    tasks.append(lv.apply(fces_app,c,csvdata))


result = [task.get() for task in tasks]  # blocks until all results are back
print('done')


################################################################################



