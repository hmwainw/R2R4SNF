###############################################################################
"""
    Last modified: August 13, 2019
    
    @author: Milos Atz <milos.atz@berkeley.edu
    
    comments: Class controlling the calculation of repository array footprint
    
"""
###############################################################################
import os
import imp
import time
import datetime
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import nwpy.repository_area
from nwpy.repository_area.data import thermalproperties
from nwpy.repository_area.heat import HeatTransfer
from nwpy.repository_area.iter import Iteration
from nwpy.repository_area.waste import Waste
###############################################################################


def fces_app(rep_inst, filedata):
    """This function is a wrapper to call the instance method calculate_app 
    in parallel calculations on Savio - it is built specifically for the
    parametric studies over the FCES analysis examples """
    
    app = rep_inst.calculate_app(log=True)
    # filedata is: [eg, str(stg), wf.form, pl, pkgs, r, st]
    txt = ','.join(filedata)+','+str(app)+'\n'
    outfile = 'fces_app_data.csv'
    out = open(os.path.join('.', outfile), 'a')
    out.writelines(txt)
    out.close()


def fces_min_st(rep_inst, filedata):
    """This function is a wrapper to call the instance method 
    calculate_required_st in parallel calculations on Savio - it is 
    built specifically for the parametric studies over the FCES analysis 
    examples """
    
    rst = rep_inst.calculate_required_st()
    # filedata is: [eg, str(stg), wf.form, pl, pkgs, r, st]
    txt = ','.join(filedata)+','+str(rst)+'\n'
    outfile = 'fces_rst_data.csv'
    out = open(os.path.join('.', outfile), 'a')
    out.writelines(txt)
    out.close()


def calculate_area_per_package(rep_inst, r, w, n, st):
    """This function is a wrapper to call the instance method
    calculate_app in parallel calculations on Savio"""
    
    app = rep_inst.calculate_app(log=True)
    text = str(st)+','+str(n)+','+str(app)+'\n'
    outfile = str(r)+'_'+str(w)+'.app'
    # structure is st, nwf, app
    out = open(os.path.join(rep_inst.wdpath, outfile), 'a')
    out.writelines(text)
    out.close()
    return(app)


def calculate_surface_storage_times(rep_inst, r, w, n):
    """Calculate the required surface storage time for possible disposal and
    for disposal with minimum repository footprint
        
    """
    
    required_st = rep_inst.calculate_required_st()
    maximum_st = rep_inst.calculate_maximum_st()
    text = str(n)+','+str(required_st)+','+str(maximum_st)+'\n'
    outfile = str(r)+'_'+str(w)+'.st'
    # structure is st, nwf, app
    out = open(os.path.join(outfile), 'a')
    out.writelines(text)
    out.close()
    return([required_st, maximum_st])


###############################################################################
# CANISTER ARRAY CLASS
###############################################################################


class Repository(object):


    """
    The Repository class represents a subsection of the repository 
    where waste of a single type is emplaced all at once. If a fuel 
    cycle produces multiple wastes, the repository will have multiple 
    sections. The Repository class instantiates Repository objects to 
    determine the area required for disposal of the waste in that 
    section and sums them to determine the total area required for the 
    fuel cycle.
    
    """


    def __init__(self, rock, w=None, st=0.0, depth=500.0,
                 surface_temp=15.0, ambient_temp=None, tol=1.0e-2,
                 test=False):
        """Set up the general attributes of a repository
        
        Parameters
        ----------
        rock: str
            granite, clay, or salt
        
        waste (optional):
            Waste to be disposed of in the repository
        
        st (optional): float
            Storage time (y) between waste discharge and emplacement
        
        depth (optional): float
            Depth (m) of repository horizon; determines ambient temp
            
        surface_temp (optional): float
            Above-ground temperature (degrees C)
            
        ambient_temp (optional): float
            Ambient temperature at repository horizon (degrees C)
            
        tol (optional): float
            Solution tolerance
            
        """
        
        # SET ATTRIBUTES
        self.name = rock+' repository array'
        self.N = 9 # array size=NxN; # adjacent drifts/packages=(N-1)/2
        self.tol = tol
        self.st = st
        self.datapath = nwpy.repository_area.get_data('')
        self.wdpath = os.path.join(os.getcwd(), rock)
        self._get_rock_thermal_properties(rock)
        if(not ambient_temp):
            self.ambient_temp=self.calc_ambient_temperature(depth,surface_temp)
        else:
            self.ambient_temp = ambient_temp
        if(w is not None):
            self.load_waste(Waste(w))
        # HEAT TRANSFER INSTANCE
        self.ht = HeatTransfer(st)
        # ITERATION INSTANCE
        self.iter = Iteration()


    def __repr__(self):
        """repr method for Repository class"""
    
        return('Repository array object')
    
    
    def __str__(self):
        """str method for Repository class"""
    
        p = self.name+' for '+self.waste.desc+' '+self.waste.id


    def _get_rock_thermal_properties(self, rocktype):
        """Import the thermal properties for relevant materials"""
        
        self.rock = {}
        self.rock['name'] = rocktype
        self.rock['k'] = thermalproperties.k[rocktype]
        self.rock['a'] = thermalproperties.a[rocktype]


    @staticmethod
    def calc_ambient_temperature(depth, above_ground_temp=15.0):
        """Calculate the ambient temperature at the repository horizon as a
            function of repository depth and above-ground temperature"""
        
        # note: above_ground_temp included as trivial optional arg
        return(above_ground_temp+25.0*depth/1000.0)


    #--------------------------------------------------------------------------
    # REPOSITORY EBS DESIGN
    # The EBS design relies on data about the host rock and the waste
    # package. Once the Repository is instantiated and the host rock is
    # specified, and once waste characteristics are specified, the design
    # of the repository can be pulled in from data files.
    #--------------------------------------------------------------------------
    def load_waste(self, waste_inst):
        """Capture all information from waste instance
        
        waste_inst: repository_area Waste instance
            Contains all relevant information for footprint calculation
            
        """
        
        self.update_wdpath(waste_inst)
        self.waste = {}
        self.waste['id'] = waste_inst.id # waste id
        self.waste['n_wf'] = waste_inst.n_wf # pkg loading
        self.decay_heat = waste_inst.decay_heat
        for x in ['evaluationgroup', 'stagenumber', 'name',
                  'form', 'number', 'loading_fraction']:
            try:
                self.waste[x] = getattr(waste_inst, x)
            except AttributeError:
                pass
        if(hasattr(waste_inst, 'batch') and 'form' in self.waste.keys()):
            self.waste['form'] += str(waste_inst.batch)
        # modify the repository design based on waste attributes
        self.design_repository(waste_inst.id)
        self.set_package(waste_inst.pkg['diameter'], waste_inst.pkg['length'])
        if(hasattr(waste_inst, 'st')):
            self.update_st(waste_inst.st)

    
    def update_st(self, st):
        """Update the surface storage time, as well as the heat transfer
            attribute that relies on it"""
        
        self.st = st
        self.ht = HeatTransfer(st)
    
    
    def update_wdpath(self, waste_inst):
        """Update the working directory path based on waste information"""
    
        if(hasattr(waste_inst, 'evaluationgroup')):
            self.wdpath = os.path.join(self.wdpath, waste_inst.evaluationgroup)
            if(hasattr(waste_inst, 'stagenumber')):
                self.wdpath=os.path.join(self.wdpath,str(waste_inst.stagenumber))
        else:
            subdir = waste_inst.id
            for x in ['name', 'form']:
                if(hasattr(waste_inst, x)):
                    subdir+='_'+getattr(waste_inst, x)
            self.wdpath = os.path.join(self.wdpath, subdir)
        
        
    def design_repository(self, wastetype):
        """Import the repository layout and design data
        
        Parameters
        ----------
        wastetype: str
            'SNF' or 'UNF' or 'HLW'; case insensitive
        
        diameter (optional): float
            package diameter (m), depends on number of waste forms
        
        length (optional) float
            package length (m)
            
        """
        
        if(wastetype.lower()=='unf'):
            wastetype = 'snf'
        self.ebs = {}
        file = os.path.join(self.datapath, self.rock['name']+'.py')
        temp = imp.load_source(self.rock['name']+'.py', file)
        for item in [a for a in dir(temp) if not a.startswith('__')]:
            if(item=='ebs'): # item is dict
                ebs_temp = getattr(temp, item)
                for key in ebs_temp[wastetype.lower()]: # item
                    self.ebs[key] = ebs_temp[wastetype.lower()][key]
            else:
                setattr(self, item, getattr(temp, item))
        self._get_ebs_thermal_properties()
    
    
    def _get_ebs_thermal_properties(self):
        """Assign thermal conductivity values to ebs layers"""
        
        self.ebs['k']=[]
        for mat in self.ebs['material']:
            self.ebs['k'].append(thermalproperties.k[mat])


    def set_package(self, diameter, length, cap_flag=False):
        """Account for package (overpack) thickness, calculate drift radius
        
        Parameters
        ----------
        diameter: float
            package diameter (m), depends on number of waste forms
            
        length: float
            package length (m)
        
        cap_flag (optional): Boolean
            Use if you want to add overpack thickness to package length
        
        """
        
        self.pkg = {}
        op_thickness = self.ebs['dr'][self.ebs['layer'].index('overpack')]
        self.pkg['d'] = diameter+2*op_thickness
        if(cap_flag): # WasteForm instance; account for op in length
            self.pkg['l'] = length+2*op_thickness
        else:
            self.pkg['l'] = length
        self.ebs['r_drift'] = sum(self.ebs['dr'][:-1])+(self.pkg['d']/2.0)


    #--------------------------------------------------------------------------
    # CALCULATION OF REQUIRED SURFACE STORAGE TIME
    # These methods are used to determine the surface storage time required
    # before disposal of a given waste type with some package loading is
    # possible in a repository. The functions search through increasing
    # values of surface storage time and utilize the "single_package" method
    # from the heat transfer transfer module, which returns a Boolean
    # indicating whether disposal is possible.
    #--------------------------------------------------------------------------
    def calculate_required_st(self, verbose=False):
        """Determine the surface storage time required before waste can be
        emplaced in the repository without breaching thermal constraints.
        
        Parameters
        ----------
        self: Repository instance
        
        verbose (optional): Boolean
            Print some stuff
        
        Returns
        -------
        Surface storage time (float) required for possible disposal
        
        """
        
        # instantiate separate heat transfer module so as not to
        # overwrite the one built into the array instance.
        st = 0.0 # starting value
        steps = [100.0, 20.0, 5.0, 1.0]
        # initial run
        ht = HeatTransfer(st)
        disp = ht.single_package(self)
        if(disp):
            return(st)
        else: # requires > 0 y storage; should be everything
            for dt in steps:
                while(disp == False):
                    st += dt
                    ht = HeatTransfer(st)
                    disp = ht.single_package(self)
                    if(verbose):
                        print('st = '+str(st)+'  temp = '+str(peak_temp))
                st = st-dt
                disp = False
            return(st+dt)


    def calculate_maximum_st(self, tol=1e-3, sensitivity=0.05, verbose=False):
        """Determine the surface storage time required before waste can
        be emplaced with the minimum possible repository dimensions.
        
        Parameters
        ----------
        self: Repository instance
        
        tol (optional): float
            Solution tolerance
            
        verbose (optional): Boolean
            Print some stuff
        
        Returns
        -------
        Surface storage time (float) required to reach minimum footprint
        
        """
        
        # instantiate separate heat transfer module so as not to
        # overwrite the one built into the array instance.
        st = 0.0 # starting value
        steps = [500, 100.0, 20.0, 5.0, 1.0]
        #error = 1.0
        # set minimum dimensions
        self.spacing['drift'] = 2*self.ebs['r_drift']
        self.spacing['pkg'] = self.pkg['l']
        for dt in steps:
            peak_temp = self.constraint['wp']+1
            err = 1.0
            while(peak_temp > self.constraint['wp']):
                st += dt
                ht = HeatTransfer(st)
                peak_temp = ht.find_peak_wp_temp(self)
                if(verbose==True):
                    print('st = '+str(st)+' N = '+str(self.N)+
                          ' temp = '+str(peak_temp))
                if(peak_temp < self.constraint['wp']): # value converged
                    while(err > sensitivity):
                        err = ht.calc_array_size_sensitivity(self, peak_temp)
                        if(verbose==True):
                            print('st = '+str(st)+' N = '+str(self.N)+' temp'+
                                  ' = '+str(peak_temp)+' err = '+str(err))
                        # increase array size, update peak temp
                        self.N += 2
                        peak_temp = peak_temp*(err+1)
            self.N -= 2 # revert array size before back stepping in time
            err = 1.0 # reset array error
            st = st - dt
        self.N = 9 # reset N
        return(st+dt)


    #--------------------------------------------------------------------------
    # REPOSITORY ARRAY FOOTPRINT CALCULATION
    # These functions carry out the iterations to determine the repository
    # footprint.
    #--------------------------------------------------------------------------
    def calculate_footprint(self):
        """Multiply the area required per package by the total
        number of packages to find the required disposal footprint
        
        Parameters
        ----------
        self: repository instance
            contains rock and waste data
        
        Returns
        -------
        Area required for disposal of all packages of the given waste (sqm)
        
        """
        
        # maybe move this outside the class and instantiate arrays inside
        area_per_pkg = self.calculate_app()
        return(area_per_pkg*np.ceil(self.waste['number']))
    
    
    def calculate_app(self, guess=None, array_sensitivity=0.05,
                      log=False, test=False):
        """Calculate minimum drift and package spacing of a square array, 
        while ensuring it is large enough to use as unit cell for a larger 
        repository.
        
        Parameters
        ----------
        self: repository instance
            contains rock and waste data
        
        guess (optional): list or tuple
            2-item list or tuple containing starting guesses for drift 
            and package spacing
            
        array_sensitivity (optional): float
            Allowable relative error for the sensitivity of the constraint 
            to packages outside the simulated array.
            
        log (optional): Boolean
            Indicates whether to set up a working directory and log outputs 
            in an output file 
            
        test (optional): Boolean
            Indicates whether the method is being tested (in which case a 
            special output directory is made
            
        Returns
        -------
        The area required per package in the repository
            
        """
        
        if(test): # always log tests
            log = True
        if(log==True and not os.path.exists(self.wdpath)):#make wd if logging
            os.makedirs(self.wdpath)
        err = 1.0
        if(not self.ht.single_package(self)):
            print('No solution; single package violates constraint')
            return(0.0)
        else:
            while(err > array_sensitivity): # find the footprint
                t0 = time.time()
                res = self.footprint_iter(guess)
                area_per_pkg = res.fun/self.N/self.N
                runtime = round((time.time()-t0)/60.0, 2) # min
                # Update repository dimensions
                self.spacing['drift'], self.spacing['pkg'] = res.x
                # test for peak rock temperature (should be met)
                #peak_rock_temp = self.ht.find_peak_rock_temp(self)
                # Test sensitivity of wp peak temp to pkg outside array
                T = self.iter.get_peak_temp(res.fun)
                err=self.ht.calc_array_size_sensitivity(self, T)
                #print err
                # Write the data
                if(log==True):
                    self.iter.write_data(self, res, err, array_sensitivity)
                # Update about results of current iteration
                self.iter.reset()
                self.N += 2
            # revert N to final value from iterations
            self.N += -2
            return(area_per_pkg)
    
    
    def footprint_iter(self, guess=None):
        """Given N, calculate minimum drift and package spacing of an NxN
        array constrained by the waste package surface temperature limit
        
        Parameters
        ----------
        self: repository instance
            contains rock and waste data
        
        guess (optional): list
            2-item list containing guesses for the drift and package spacing
        
        """
        
        # constraints
        cons = ({'type':'ineq', 'fun':lambda x: (self.constraint['wp']-\
                                                 self.calc_peak_temp(x))},
                {'type':'ineq', 'fun':lambda x: x[0] - 2*self.ebs['r_drift']},
                {'type':'ineq', 'fun':lambda x: x[1] - self.pkg['l']},
                {'type':'ineq', 'fun':lambda x: x[0]-x[1]})
        # arguments
        fxn_args = [self.ebs['r_drift'], self.pkg['l'], self.N]
        # assign guess depending on what's already been calculated
        if(self.iter.array_idx==0):
            if(guess):
                g = guess
            else:
                g = [min(self.spacing.values()), min(self.spacing.values())]
        else:
            g=self.spacing.values()
        # perform the minimization
        return(minimize(self._calc_area, g, args=fxn_args,
                        method='COBYLA', constraints=cons, tol=self.tol))
    
    
    def calc_peak_temp(self, dims=None, iter=True):
        """Based on drift and package spacing, calculates waste package peak 
        temperature at center of canister array
        
        Parameters
        ----------
        self: repository instance
            contains rock and waste data
    
        spacing (optional): list
            2-item list containing new drift and package spacing
    
        """
        
        if(dims is not None): # update drift, package spacing
            self.spacing['drift'], self.spacing['pkg'] = dims
        T = self.ht.find_peak_wp_temp(self)
        A = self._calc_area([self.spacing['drift'], self.spacing['pkg']],
                            [self.ebs['r_drift'],self.pkg['l'],self.N])
        if(iter):
            self.iter.read(ds=self.spacing['drift'], wps=self.spacing['pkg'],
                           area=A, temp=T)
            self.iter.update()
            print(self.iter)
        return(T)


    @staticmethod
    def _calc_area(spacing, args):
        """Calculate the footprint of the repository given its dimensions
            
        Parameters
        ----------
        spacing: list
            Contains the repository drift and package spacing
        
        args: list
            1. drift radius
            2. the number of drifts in the array
            3. package length
            4. the number of packages per drift
        
        Returns
        -------
        Repository footprint (float; units: square meters)
        
        """
        
        drift_spc, pkg_spc = spacing
        drift_r, pkg_len, N = args
        # Repository length is parallel to drift; must account for extra
        # half-package length on the ends for last packages in drift
        l = (N-1)*pkg_spc
        l += 2*(0.5*pkg_len)
        # Repository width is perpendicular to drift; account for drift
        # radius on the outsides of the first and last drifts
        w = (N-1)*drift_spc
        w += 2*(drift_r)
        return(l*w)


###############################################################################

