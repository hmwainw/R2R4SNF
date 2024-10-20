###############################################################################
"""
    Last edited on August 13, 2019

    @author: matz

    comments: Load required waste information into waste instance
    
"""
###############################################################################
import os
import copy
import numpy as np
import pandas as pd
import scipy
from scipy.interpolate import interp1d
###############################################################################


class Waste(object):


    def __init__(self, data_obj, st=None):
        """The properties of the waste are loaded into a Waste instance for
            the footprint calculation. These can be loaded either from a 
            WasteForm instance or from a path to a CSV file supplied by the user.
            
        The required information is as follows:
        - id: SNF or HLW
        - n_wf: Number of waste forms
        - package diameter
        - package length
        - decay heat characteristics; heat (W) tabulated as fxn of time (y)
        
        Optional information:
        - Evaluation group
        - Stage number
        - Waste name (e.g. "uox")
        - Waste form (e.g. "metal" vs. "ceramic")
        
        Parameters
        ----------
        dataobj: path to datafile OR nwpy WasteForm instance
        
        """
        
        if(st is not None):
            self.st = st
        if(type(data_obj)==str): # treat data_obj as path
            self.load_from_file(data_obj)
        else: # treat data_obj as WasteForm instance
            self.load_from_stream(data_obj)


    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def load_from_stream(self, wf_inst):
        """Load waste from WasteForm instance
            
        Parameters
        ----------
        self: Waste instance
        
        wf_inst: WasteForm instance
            Contains information added as attributes to Repository instance
        
        st (optional): float
            Surface storage time for waste
            
        """
        
        w = copy.deepcopy(wf_inst) # make temp obj
        time = [float(x) for x in list(w.heat.columns)]
        heat = [sum(w.heat[t]) for t in w.heat]
        # check for optional attributes to pass
        opt_attr = {}
        for attr in ['evaluationgroup','stagenumber','name',
                     'form','number','batch', 'loading_fraction']:
            try:
                opt_attr[attr] = getattr(w, attr)
            except AttributeError: # waste doesn't have this attribute dummy
                pass
        if(not hasattr(w, 'loading')):
            w.loading = 1 # hlw where it may not be specified
        self.load_waste(w.id, w.loading, w.canister['Diameter'],
                        w.canister['Length'], time, heat, **opt_attr)

    
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    def load_from_file(self, data_file_path):
        """Load the decay heat and package data from a csv file
            
        Parameters
        ----------
        self: Waste instance
        
        data_file_path: str
            Path to csv containing waste data
            
        st (optional): float
            Surface storage time for waste
            
        """
        # read heat and kwarg data from table
        x = pd.read_csv(data_file_path, index_col=0, header=None)
        time, heat = [], []
        for key in x.index:
            try:
                time.append(float(key))
                heat.append(float(x.loc[key][1]))
            except ValueError: # cannot convert str to float: not heat data
                pass
        # pull optional data from x
        opt_args = {}
        for attr in ['evaluationgroup','stagenumber','name','form','number']:
            try:
                opt_args[attr] = x.loc[attr][1]
            except KeyError: # datafile doesn't have this attribute dummy
                pass
        # number = float(x.loc['number'][1])
        self.load_waste(x.loc['id'][1], float(x.loc['n_wf'][1]),
                        float(x.loc['pkg_diameter'][1]),
                        float(x.loc['pkg_length'][1]),
                        time, heat, **opt_args)


    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def load_waste(self, id, n_wf, pkg_d, pkg_l, time, heat, **kwargs):
        """Load the information taken from the WasteForm instance or the data
        file into the Repository instance

        Parameters
        ----------
        self: Repository instance

        id: str
            Waste identifier (snf vs. hlw)

        n_wf: int
            Number of waste forms per package

        pkg_d: float
            Diameter of waste package (m)

        pkg_l: float
            Length of waste package (m)

        time: Sequence (list or tuple)
            Data points for decay heat interpolation (time, y)

        heat: Sequence (list or tuple)
            Data points for decay heat interpolation (heat, W)

        kwargs
        ------
        - name: str
            e.g. "uox"
        - form: str
            e.g. "metal" or "ceramic"
        - evaluationgroup: str
            e.g. "eg01"
        - stagenumber: str
            e.g. "1" corresponds to EGXX-1

        """

        # assign data to waste dictionary
        self.id = id
        self.n_wf = n_wf
        for x in kwargs.keys(): # run through optional args
            setattr(self, x, kwargs[x])
        self.qinterp = self._make_heat_interpolation(time, heat)
        self.pkg = {'diameter': pkg_d, 'length': pkg_l}


    #--------------------------------------------------------------------------
    # DECAY HEAT
    # Based on an interpolation of decay heat data, this function returns
    # the decay heat (W) generated from a waste package. This heat function
    # is used as the source term in the heat conduction calculations to
    # evaluate the repository footprint and temperature constraints.
    #--------------------------------------------------------------------------
    def decay_heat(self, time):
        """Return the thermal power of the waste package at a given time"""
            
        if(not time > 0):
            raise ValueError('Time must be greater than 0')
        return(10.0**self.qinterp(time))


    @staticmethod
    def _make_heat_interpolation(time, heat):
        """Interpolate the heat data when waste is loaded"""
            
        logq = np.zeros(len(heat))
        for i in range(0, len(heat)):
            logq[i] = np.log10(heat[i])
        return(interp1d(time, logq, fill_value="extrapolate"))


###############################################################################
