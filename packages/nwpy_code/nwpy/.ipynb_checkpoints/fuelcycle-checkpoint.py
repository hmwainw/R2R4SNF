###################################################################################
"""
    Last edited on June 28, 2019

    @author: matz

    comments: The FuelCycle class instantiates and executes Stage instances to 
              produce waste streams for all stages in a fuel cycle
    
"""
###################################################################################
import os
import imp
import decimal
import numpy as np
import pandas as pd
import scipy
from scipy import optimize
import matplotlib.pyplot as plt
from collections import Counter
import nwpy
from nwpy import stage
###################################################################################


class FuelCycle(object):


    """
    A fuel cycle is a collection of stages, which are each comprised of an
    irradiation system and the processes that support it, such as fuel
    fabrication and used fuel reprocessing. The FuelCycle object is used to
    perform calculations on all stages of an evaluation group at once. In
    particular, it can be used to (a) quickly generate data for the main 
    wastes produced in that fuel cycle, and (b) benchmark that data against
    what was reported in the DOE Fuel Cycle Evaluation and Screening study
    Metric Data (Appendix D).
    
    """


    def __init__(self, evaluationgroup, tol=1e-9, test=False):
        """Initialize the FuelCycle class"""
        
        # SET ATTRIBUTES
        self.name =  evaluationgroup
        self.evaluationgroup = evaluationgroup
        self.test = test
        # DIRECTORY SETUP
        # absolute path to data directory
        maindir = os.getcwd()
        #print(os.path.abspath('./'))
        # datapath
        #self.datapath = get_data('') #maindir+'/data/'
        self.datapath = nwpy.get_data('', test) #maindir+'/data/'
        # set up output dir: first, test that output dir is made
        outdir = maindir+'/output/'
        if(not os.path.isdir(outdir)):
            os.mkdir(outdir)
        if(test): # IF TESTING, WORKING DIR WILL BE IN output/test/
            outdir += 'test/'
            if(not os.path.isdir(outdir)):
                os.mkdir(outdir)
        # then, test that evaluation group dir is made
        outdir += self.evaluationgroup
        if(not os.path.isdir(outdir)):
            os.mkdir(outdir)
        # find out number of stages; load the fuel cycle data
        self._get_total_stages()
    
    
    def __repr__(self):
        return('FuelCycle instance: '+str(self.name))
    
    
    def __str__(self):
        fc_datapath = self.datapath+'fc/'+str(self.evaluationgroup)+'.fc'
        try:
            fcd = open(fc_datapath).read().splitlines()
        except IOError as error:
            print('Fuel cycle data file for '+self.name+\
                  ' does not exist: '+fc_datapath)
            raise
        p = ''
        for l in fcd:
            line = l.split()
            if(int(line[0])==0):
                if(line[1]=='Item'):
                    p += 'STG'+l[1:]+'\n'
                else:
                    p += ' '.join(line[1:])+'\n'
            else:
                p += ' '.join(l[0:2])+' '+l[2:]+'\n'
        return(p[:-1])

        
    def _get_total_stages(self):
        """Find number of stages in the fuel cycle"""
        
        fc_datapath = self.datapath+'fc/'+str(self.evaluationgroup)+'.fc'
        try:
            fcd = open(fc_datapath).read().splitlines()
        except IOError as error:
            print('Fuel cycle data file for '+self.name+\
                  ' does not exist: '+fc_datapath)
            raise
        self.totalstages = max([int(line[0]) for line in fcd])
    
    
    #------------------------------------------------------------------------------
    # BENCHMARKING
    # To confirm that the fuel cycles are modeled correctly and that the data are
    # treated properly, we can compare some broad numerical values to data from
    # the DOE FCES study. For each evaluation group, that study produced waste
    # management data that was used to compare the fuel cycles. The relevant data
    # to this analysis is:
    # 1. Mass of SNF+HLW
    # 2. Activity of SNF+HLW at 100 years
    # 3. Activity of SNF+HLW at 100,000 years
    # For each fuel cycle, these properties are calculated and the ratio between
    # the results and the data from the DOE FCES study is returned.
    #------------------------------------------------------------------------------
    def benchmark(self, **kwargs):
        """Loop over all stages in a fuel cycle and compare SNF+HLW mass and 
        activity to the results published in the DOE FCES study
        
        Parameters
        ----------
        kwargs: Any keyword arguments to pass to accepting methods within the
        stage module (e.g. reprocess, etc).
        
        Returns
        -------
        None, but prints the ratio between calculated values and FCES metric data
        
        """
        
        fces_data = pd.read_csv(self.datapath+'metric_data.csv',index_col=0)
        fces_data = fces_data.loc[self.evaluationgroup.upper()]
        wastemass = 0.0
        activity = Counter({1e2: 0.0, 1e5: 0.0})
        for stg in range(1, self.totalstages+1):
            s = stage.Stage(self.evaluationgroup, stg, test=self.test)
            m_stg, a_stg = s.benchmark_stage(**kwargs)
            wastemass += m_stg
            activity += a_stg
        # mass case
        wastemass = wastemass/1e6/100 # g -> t/GWe-y
        fces_mass = fces_data['Mass SNF+HLW (t/GWe-y)']
        ratio = fces_data['Mass Renormalization Factor']
        print('Mass Ratio: '+str(wastemass*ratio/fces_mass))
        # radioactivity cases
        # need to figure out how to deal with stream objects AND lists of stuff
        a100 = fces_data['Activity of SNF+HLW at 100y (Ci/GWe-y)']
        a1e5 = fces_data['Activity of SNF+HLW at 1e5y (Ci/GWe-y)']
        print('Activity Ratio (100 y): '+str(activity[1e2]*ratio/a100/100.0))
        print('Activity Ratio (1e5 y): '+str(activity[1e5]*ratio/a1e5/100.0))


    #------------------------------------------------------------------------------
    # FUEL CYCLE WASTES
    # Discharge all the wastes from each fuel cycle stage, producing an array of
    # WasteForm objects that can be used by other programs.
    #------------------------------------------------------------------------------
    def discharge_wastes(self, endtime=1e5, steps=None, **kwargs):
        """Loop over all stages in a fuel cycle and return waste forms for each.
            
        Parameters
        ----------
        endtime (optional): float
            time at the end of the decay calculation time range
        
        steps (optional): int
            number of steps required for the calculation
        
        kwargs: waste loading keyword arguments
        - verbose: print information about loading
        - loading: SNF loading into packages
        - loading_fraction: HLW loading into waste forms
        
        #recycle: if applicable, recycle salt to concentrate the waste.
        #plot: if applicable, produce a plot of the loading constraints.
        #loading_level: if applicable, 'hi' or 'low' for htgr snf.
        
        Results
        -------
        Dict with stage numbers as keys containing waste form objects produced 
        by each stage
        
        """
        
        w_dict = {}
        for stg in range(1, self.totalstages+1):
            s = stage.Stage(self.evaluationgroup, stg, test=self.test)
            w_dict[stg] = s.discharge_all_wastes(endtime, steps, **kwargs)
        return(w_dict)


###################################################################################