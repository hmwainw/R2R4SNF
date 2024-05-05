###############################################################################
"""
    Last modified: July 9, 2019
    
    @author: Milos Atz <milos.atz@berkeley.edu>
    
    comments: Simple test case for error demonstration
    
"""
###############################################################################
import types
import copy_reg
import ipyparallel as ipp
import pandas as pd
import nwpy
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
# REPOSITORY TYPE
# The repository host rock and disposal concept informs the heat transfer
# properties of the rock and engineered barrier system materials, as well as
# the thermal constraint that must be met.
#------------------------------------------------------------------------------
rocktypes = ['granite', 'clay', 'salt']


#------------------------------------------------------------------------------
# Example waste data
# Read composition and heat data from files
comp = pd.read_csv('./ex_comp.csv', index_col=0)
comp.columns = pd.to_numeric(comp.columns)
heat = pd.read_csv('./ex_heat.csv', index_col=0)
heat.columns = pd.to_numeric(heat.columns)
# Assign data attributes to waste stream
mass = 1839720.5588822358
canister = {'Diameter': 0.82, 'Length': 5.0, 'Thickness': 0.125}
number = 1002
kwdata = {'loading': 4,
          'evaluationgroup': 'example-01',
          'stagenumber': 1,
          'id': 'snf',
          'form': 'snf',
          'heat': heat,
          'comp': comp}


wf = nwpy.stream.WasteForm(mass, number, canister, **kwdata)


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


###############################################################################


client = ipp.Client()
lv = client.load_balanced_view()
tasks = []


for r in rocktypes:
    c = repository_area.repository.Repository(r, wf)
    # I/O DATA: EG, STG, WF, LOAD, N_PKGS, ROCK, ST
    csvdata = [wf.evaluationgroup, str(wf.stagenumber), wf.form,
               str(wf.loading), str(wf.number), r]
    tasks.append(lv.apply(fces_min_st,c,csvdata))
    for st in storage_time:
        c.update_st(st)
        tasks.append(lv.apply(fces_app,c,csvdata+[str(st)]))


result = [task.get() for task in tasks]  # blocks until all results are back
print('done')


################################################################################



