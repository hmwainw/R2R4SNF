###############################################################################
"""
    Last edited on July 3, 2019
    
    @author: matz
    
    comments: Run file for parallel parametric FOM calculations 
        for the FCES fuel cycles
        
"""
###############################################################################
import os
import types
import copy_reg
import ipyparallel as ipp
from collections import Sequence
import nwpy
from nwpy import attractiveness
from nwpy.attractiveness import materials
from nwpy.attractiveness.criticalmass import CriticalMass
from nwpy.attractiveness.criticalmass import run as run_cm
from nwpy.attractiveness.heat import Heat
from nwpy.attractiveness.heat import run as run_heat
from nwpy.attractiveness.dose import Dose
from nwpy.attractiveness.dose import run as run_dose
###############################################################################


fuelcycles = ['eg04','eg05','eg06','eg07','eg08','eg09','eg10',
              'eg11','eg12','eg13','eg14','eg15','eg16','eg17','eg18','eg19','eg20',
              'eg21','eg22','eg23','eg24','eg25','eg26','eg27','eg28','eg29','eg30',
              'eg31','eg32','eg33','eg34','eg35','eg36','eg37','eg38','eg39','eg40']
targets = ['pu', 'tru']


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


###############################################################################


#------------------------------------------------------------------------------
# Parallelization
#------------------------------------------------------------------------------
client = ipp.Client()
lv = client.load_balanced_view()


#------------------------------------------------------------------------------
# GENERATE WASTE STREAMS - store in dictionary
#------------------------------------------------------------------------------
fcw = {}
for eg in fuelcycles:
    fcw[eg] = {}
    fc = nwpy.fuelcycle.FuelCycle(eg)
    for stg in xrange(1, fc.totalstages+1):
        fcw[eg][stg] = [] # append waste streams to this list
        s = nwpy.stage.Stage(eg, stg)
        strm = s.discharge_streams()
        strm = s.cool(strm)
        strm = s.reprocess(strm)
        wf = s.load_waste(strm) # use default loading
        if(not isinstance(wf, Sequence)):
            wf = [wf]
        # Modify waste stream to account for inert media
        for wfi in wf:
            if(wfi.mass_fraction('actinide')==0.0): # skip
                continue
            else:
                wfi = materials.merge(s, wfi)
                wfi = s.decay(wfi, endtime=1e6, steps=10)
                fcw[eg][stg].append(wfi)


print('Waste stream generation: done')


#------------------------------------------------------------------------------
# DOSE CALCULATION - FIRST PARALLELIZATION
#------------------------------------------------------------------------------
tasks = []
for eg in fcw:
    for stg in fcw[eg]:
        for wfi in fcw[eg][stg]:
            print(eg+'-'+str(stg)+wfi.form)
            idx = 0
            for time in wfi.comp.columns:
                if(time < 5.0):
                    continue
                else:
                    wdpath = os.path.join('.','output', eg, str(stg))
                    # DOSE ----------------------------------------------------
                    wd_dose = os.path.join(wdpath, 'dose','t'+str(idx))
                    a = [wfi, time, wd_dose, idx]
                    tasks.append(lv.apply(run_dose, Dose, a))
                # increase time index
                idx+=1


result = [task.get() for task in tasks]  # blocks until all results are back
print('Dose results: done')


#------------------------------------------------------------------------------
# CRITICAL MASS AND HEAT CALCULATIONS - SECOND PARALLELIZATION
#------------------------------------------------------------------------------
tasks = []
for eg in fcw:
    for stg in fcw[eg]:
        for wfi in fcw[eg][stg]:
            print(eg+'-'+str(stg)+wfi.form)
            idx = 0
            for time in wfi.comp.columns:
                if(time < 5.0):
                    continue
                else:
                    wdpath = os.path.join('.','output', eg, str(stg))
                    for target in targets:
                        # Critical mass
                        wd_cm = os.path.join(wdpath, 'cm','t'+str(idx))
                        b = [wfi, time, target, wd_cm, idx]
                        tasks.append(lv.apply(run_cm, CriticalMass, b))
                        # Heat
                        wd_heat = os.path.join(wdpath, 'heat','t'+str(idx))
                        c = [wfi, time, target, wd_heat, idx]
                        tasks.append(lv.apply(run_heat, Heat, c))
                # increase time index
                idx+=1


result = [task.get() for task in tasks]  # blocks until all results are back
print('Heat, critical mass results: done')


###############################################################################

