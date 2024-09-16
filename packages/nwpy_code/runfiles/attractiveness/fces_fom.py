###############################################################################
"""
    Last edited on July 16, 2019
    
    @author: matz
    
    comments: Run file for parallel parametric FOM calculations 
        for the FCES fuel cycles
        
"""
###############################################################################
import os
import time
import types
import pickle
import copy_reg
import ipyparallel as ipp
from collections import Sequence
import nwpy
from nwpy import attractiveness
from nwpy.attractiveness import materials
from nwpy.attractiveness.criticalmass import CriticalMass
from nwpy.attractiveness.criticalmass import run as run_cm
#from nwpy.attractiveness.heat import Heat
#from nwpy.attractiveness.heat import run as run_heat
from nwpy.attractiveness.dose import Dose
from nwpy.attractiveness.dose import run as run_dose
###############################################################################
fuelcycles = ['eg01','eg02','eg03','eg04','eg05','eg06','eg07','eg08','eg09',
              'eg10','eg11','eg12','eg13','eg14','eg15','eg16','eg17','eg18',
              'eg19','eg20','eg21','eg22','eg23','eg24','eg25','eg26','eg27',
              'eg28','eg29','eg30','eg31','eg32','eg33','eg34','eg35','eg36',
              'eg37','eg38','eg39','eg40']
targets = ['u', 'pu', 'tru']


def get_heat(w_inst, target, time):
    h=0.0
    for iso in w_inst.heat[time].index:
        if(nwpy.nuclidedata.is_group(iso, target)):
            h+= w_inst.heat[time][iso]
    h = h/(w_inst.mass*w_inst.mass_fraction(target)/1e3)
    dat = [w_inst.evaluationgroup, str(w_inst.stagenumber), w_inst.form,
           str(w_inst.mass/1e3),str(w_inst.number), str(time), target, str(h)]
    txt = ','.join(dat)+'\n'
    outfile = 'fces_fom_heat.csv'
    out = open(os.path.join('.', outfile), 'a')
    out.writelines(txt)
    out.close()


def get_wf_props(wf, time):
    dat = [wf.evaluationgroup, str(wf.stagenumber), wf.form,
           str(wf.mass/1e3),str(wf.number), str(time)]
    # PLUTONIUM
    pu = wf.mass_fraction('pu')*wf.mass # g/pkg
    fpu_iso = ['pu239', 'pu241']
    fpu = sum([wf.mass_fraction(x) for x in fpu_iso])/wf.mass_fraction('pu')
    # TRANSURANICS
    tru = wf.mass_fraction('tru')*wf.mass # g/pkg
    ftru_iso = ['pu239', 'pu241', 'am242m', 'cm243', 'cm245', 'cm247']
    ftru=sum([wf.mass_fraction(x) for x in ftru_iso])/wf.mass_fraction('tru')
    # URANIUM
    u = wf.mass_fraction('u')*wf.mass # g/pkg
    fu_iso = ['u233', 'u235']
    fu = sum([wf.mass_fraction(x) for x in fu_iso])/wf.mass_fraction('u')
    fu233 = wf.mass_fraction('u233')/wf.mass_fraction('u')
    dat = dat+[str(u),str(fu),str(fu233),str(pu),str(fpu),str(tru),str(ftru)]
    txt = ','.join(dat)+'\n'
    outfile = 'fces_strm_props.csv'
    out = open(os.path.join('.', outfile), 'a')
    out.writelines(txt)
    out.close()


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
# Methods to save and load nwpy.stream objects
#------------------------------------------------------------------------------
def make_case_name(strm):
    stg_name = strm.evaluationgroup+'-'+str(strm.stagenumber)
    name = nwpy.origen.Origen.make_case(strm, stg_name)
    return(name)


def save_strm(strm):
    name = make_case_name(strm)
    with open('strms/'+name+'.pkl', 'wb') as f:
        pickle.dump(strm, f, pickle.HIGHEST_PROTOCOL)


def load_strm(name):
    with open('strms/'+name+'.pkl', 'rb') as f:
        return pickle.load(f)


###############################################################################
# Parallelization
#------------------------------------------------------------------------------
client = ipp.Client()
lv = client.load_balanced_view()
tasks = []


#------------------------------------------------------------------------------
# GENERATE WASTE STREAMS - store in dictionary
#------------------------------------------------------------------------------
fcw = {}
os.mkdir('./strms')
for eg in fuelcycles:
    fcw[eg] = {}
    fc = nwpy.fuelcycle.FuelCycle(eg)
    for stg in xrange(1, fc.totalstages+1):
        fcw[eg][stg] = [] # append waste streams to this list
        s = nwpy.stage.Stage(eg, stg)
        strm = s.discharge_streams()
        strm = s.cool(strm)
        strm = s.reprocess(strm)
        if(eg in ['eg15', 'eg16'] and stg==2):
            wf = s.load_waste(strm, loading=3)
        else: # use default loading
            wf = s.load_waste(strm)
        if(not isinstance(wf, Sequence)):
            wf = [wf]
        for wfi in wf: # skip if no actinides
            if(wfi.mass_fraction('actinide')==0.0): # skip
                continue
            else:
                wfi = materials.merge(s, wfi) # account for inert media
                wfi = s.decay(wfi, endtime=1e6, steps=10)
                fcw[eg][stg].append(wfi)
                # and for good measure, save them
                save_strm(wfi)
                #--------------------------------------------------------------
                # HEAT
                #--------------------------------------------------------------
                for t in wfi.comp.columns:
                    if(time < 5.0):
                        continue
                    else:
                        get_wf_props(wfi, t)
                        for target in targets:
                            get_heat(wfi, target, t)


print('Waste stream generation: done')
print('Heat content: done')


#------------------------------------------------------------------------------
# CRITICAL MASS - FIRST PARALLELIZATION
#------------------------------------------------------------------------------
tasks = []
for eg in fcw:
    for stg in fcw[eg]:
        for wfi in fcw[eg][stg]:
            print(eg+'-'+str(stg)+wfi.form)
            idx = 0
            for t in wfi.comp.columns:
                if(t < 5.0):
                    continue
                else:
                    wdpath = os.path.join('.','output', eg, str(stg))
                    for target in targets:
                        # Critical mass
                        wd_cm = os.path.join(wdpath, 'cm','t'+str(idx))
                        name=nwpy.origen.Origen.make_case(wfi,eg+'-'+str(stg))
                        name += '_t'+str(idx)+'_'+target
                        b = [wfi, t, target, wd_cm, name]
                        tasks.append(lv.apply(run_cm, CriticalMass, b))
                # increase time index
                idx+=1


result = [task.get() for task in tasks]  # blocks until all results are back
print('Critical mass results: done')


#------------------------------------------------------------------------------
# DOSE - SECOND PARALLELIZATION
#------------------------------------------------------------------------------
tasks = []
for eg in fcw:
    for stg in fcw[eg]:
        for wfi in fcw[eg][stg]:
            name=nwpy.origen.Origen.make_case(wfi,eg+'-'+str(stg))
            print(eg+'-'+str(stg)+wfi.form)
            idx = 0
            for t in wfi.comp.columns:
                if(t < 5.0):
                    continue
                else:
                    wdpath = os.path.join('.','output', eg, str(stg))
                    wd_dose = os.path.join(wdpath, 'dose','t'+str(idx))
                    a = [wfi, t, wd_dose, name+'_t'+str(idx)]
                    tasks.append(lv.apply(run_dose, Dose, a))
                    time.sleep(15)
                # increase time index
                idx+=1


result = [task.get() for task in tasks]  # blocks until all results are back
print('Dose results: done')


###############################################################################

