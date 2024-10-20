###############################################################################
"""
    Last edited on June 28, 2019

    @author: matz

    comments: The Stage class contains the methods to study material streams
              in the back end of the DOE FCES fuel cycle analysis examples
    
"""
###############################################################################
import os
import numpy as np
import pandas as pd
from collections import Counter
from collections.abc import Sequence
# import __init__
import nwpy
from nwpy import origen
from nwpy import stream
from nwpy import reactor
from nwpy import separation
from nwpy import loading
###############################################################################


class Stage(object):


    """
    In the context of the fuel cycle, a stage is defined as an 
    irradiation system and its supporting infrastructure, including 
    fuel manufacturing and processing operations. The Stage class is 
    used to hold the data about the stage and its associated material 
    streams, as well as the methods that act on those streams. A stage 
    instance is defined within the context of fuel cycle evaluation 
    group.
    
    The data held in the stage are required for the determination of 
    the mass and composition of wastes. Almost all data is taken from 
    the Department of Energy (DOE) Nuclear Fuel Cycle Evaluation and 
    Screening Study (FCES; 2014). Where that work does not supply data, 
    it is generated based on literature review and assumption.
    
    The methods in the Stage class reflect the actions of stage 
    operations on material streams:
    - discharge (from the irradiation sytem), 
    - cooling (after discharge, using decay cases in ORIGEN)
    - separations (if applicable), and 
    - loading of waste into disposal canisters. 
    
    These methods act on and output Stream objects, which represent 
    material streams throughout the fuel cycle. Stream objects contain 
    the data for mass and composition, as well as the decay heat 
    generation.
    
    This code should be used to generate ORIGEN decay-case input files 
    for wasteforms associated with the fuel cycles described in the DOE 
    FCES, which can be used to perform other nuclear waste management 
    calculations to compare the performance of wastes produced in 
    different fuel cycles.
    
    """
    
    
    def __init__(self, evaluationgroup, number, tol=1e-9, test=False):
        """Initialize an instance of the Stage class based on its FuelCycle"""
        
        # SET ATTRIBUTES
        self.name = evaluationgroup+'-'+str(number)
        self.evaluationgroup = evaluationgroup
        self.number = number
        self.tol = tol
        # DIRECTORY SETUP
        # absolute path to data directory
        maindir = os.getcwd()
        # datapath
        self.datapath = nwpy.get_data('', test) #maindir+'/data/'
        # set up output dir: first, test that output dir is made
        ### FIX THESE DIR SEPS
        outdir = os.path.join(maindir,'output')
        if(not os.path.isdir(outdir)):
            os.mkdir(outdir)
        if(test): # IF TESTING, WORKING DIR WILL BE IN output/test/
            outdir = os.path.join(outdir, 'test')
            if(not os.path.isdir(outdir)):
                os.mkdir(outdir)
        # then, test that evaluation group dir is made
        outdir = os.path.join(outdir, self.evaluationgroup)
        if(not os.path.isdir(outdir)):
            os.mkdir(outdir)
        # set up wd within eg dir
        self.wdpath = os.path.join(outdir, str(self.number))
        if(not os.path.isdir(self.wdpath)):
            os.mkdir(self.wdpath)
        # find out number of stages; load the fuel cycle data
        fc_datapath = os.path.join(self.datapath,'fc',
                                   str(self.evaluationgroup)+'.fc')
        self._get_stage_data(fc_datapath)
        # Initialize reactor
        self.reactor = self._get_reactor()
        # Initialize separation
        if(self.data['reprocessing']!='none'):
            self.separation = separation.Separation(datapath=self.datapath,
                                                    data=self.data,
                                                    tol=self.tol)
        # Initialize loading
        # TBD

    
    def __repr__(self):
        """repr method for the Stage class"""
        
        return('Stage instance: '+str(self.name))
    
    
    def __str__(self):
        """str method for the Stage class"""
        
        p = 'Stage instance: '+str(self.name)+'\n'
        for key in self.data:
            p = p+key+' '+str(self.data[key])+'\n'
        return(p[:-1])
    
    
    def _get_stage_data(self, fc_datapath):
        """Parse the fuel cycle data to obtain data for only the 
        selected stage"""
        
        try:
            fcd = open(fc_datapath).read().splitlines()
        except IOError as error:
            print('Fuel cycle data file for '+self.evaluationgroup+\
                  ' does not exist at path: '+fc_datapath)
            raise
        fcd = [line.split() for line in fcd]
        stgdata = [l[1:] for l in fcd if int(l[0])==self.number]
        if(stgdata == []):
            raise ValueError('Fuel cycle '+self.evaluationgroup+' has no '+
                             'stage '+ str(self.number))
        self.data = {}
        for l in stgdata:
            try: # entry is a single float
                self.data[l[0].lower()]=float(l[1])
            except:
                try: # list of strings, separated by commas
                    if('.' in l[1]): # entries are floats
                        self.data[l[0].lower()]=[float(x) for x in
                                                 l[1].split(',')]
                    else: # entries are ints
                        self.data[l[0].lower()]=[int(x) for x in
                                                 l[1].split(',')]
                except: # don't know what it is, just assign it as is
                    self.data[l[0].lower()]=l[1]


    def _get_reactor(self):
        """Factory function to instantiate reactor object; the Reactor class
        requires the stage data, the path to the data directory, and the 
        evaluation group and stage number."""

        if(self.data['reactor']=='msr'):
            return(reactor.MSR(data=self.data, datapath=self.datapath,
                               evaluationgroup=self.evaluationgroup,
                               number=self.number))
        elif(self.data['reactor']=='ffh' and
             'salttreatment' in self.data.keys()):
            return(reactor.FFH(data=self.data, datapath=self.datapath,
                               evaluationgroup=self.evaluationgroup,
                               number=self.number))
        elif(self.data['reactor']=='sfr' and 'batches' in self.data.keys()):
            return(reactor.BnBSFR(data=self.data, datapath=self.datapath,
                                  evaluationgroup=self.evaluationgroup,
                                  number=self.number))
        elif('driver' and 'blanket' in self.data.keys()):
            return(reactor.DFBB(data=self.data, datapath=self.datapath,
                                evaluationgroup=self.evaluationgroup,
                                number=self.number))
        else: # general case
            return(reactor.Reactor(data=self.data, datapath=self.datapath,
                                   evaluationgroup=self.evaluationgroup,
                                   number=self.number))


    def _get_separation(self):
        """Factory function to instantiate separation object; the Separation 
        class requires the stage data, the path to the data directory, and the
        stage mass cutoff tolerance."""
        
        if('salttreatment' in self.data.keys()):
            return(separation.LiquidFuelSep(datapath=self.datapath,
                                            data=self.data, tol=self.tol))
        else:
            return(separation.SolidFuelSep(datapath=self.datapath,
                                           data=self.data, tol=self.tol))
    
    
    #--------------------------------------------------------------------------
    # FUEL CYCLE STAGE WASTES
    # Discharge all the wastes in the fuel cycle stage, producing an array of
    # WasteForm objects that can be used by other programs.
    #--------------------------------------------------------------------------
    def discharge_all_wastes(self, endtime=1e5, steps=None, **kwargs):
        """Discharge a list of all waste streams in the stage
        
        Parameters
        ----------
        self: Stage instance
        
        endtime (optional): float
            time at the end of the decay calculation time range
            
        steps (optional): int
            number of steps required for the calculation
        
        kwargs: waste loading keyword arguments
            verbose: print information about loading .
            recycle: if applicable, recycle salt to concentrate the waste.
            plot: if applicable, produce a plot of the loading constraints.
            loading_level: if applicable, 'hi' or 'low' for htgr snf.
            
        Results
        -------
        List containing waste form objects produced by the stage
        
        """
        
        strm = self.discharge_streams()
        strm = self.cool(strm)
        if(self.data['reprocessing']!='none'):
            waste = self.reprocess(strm)
        else:
            waste = strm
        loaded_wastes = []
        if(not isinstance(waste, Sequence)):
            waste = [waste]
        for w in waste:
            try:
                temp = self.load_waste(w, **kwargs) # loading level, etc
            except:
                continue
            else:
                if(steps==None):
                    steps = np.ceil(np.log10(endtime)*20)
                temp = self.decay(temp, endtime=endtime, steps=steps)
                loaded_wastes.append(temp)
        if(len(loaded_wastes)==1):
            return(loaded_wastes[0])
        else:
            return(loaded_wastes)

    
    #--------------------------------------------------------------------------
    # BENCHMARKING
    # For each evaluation group, the DOE FCES study produced waste data that
    # to compare fuel cycles. The relevant data for this analysis is (1) Mass
    # of SNF+HLW; (2) Activity of SNF+HLW at 100 years; (3) Activity of
    # SNF+HLW at 100,000 years. For each fuel cycle stage, these properties
    # are calculated. In the benchmark function in the FuelCycle class, the
    # data from each stage are summed and the ratio between the results and
    # the data from the DOE FCES study is returned.
    #--------------------------------------------------------------------------
    def benchmark_stage(self, **kwargs):
        """Calculate the mass and activity of SNF+HLW (up through 
        reprocessing) in stage to compare with DOE FCES metric data.
        
        Parameters
        ----------
        kwargs: any argument that can be passed to any of the submethods
        called in this procedure (e.g. reprocessing, etc).
        
        
        
        Returns
        -------
        1. Mass of wastes (SNF+HLW) from stage
        2. Dict of activity of SNF+HLW at 100 and 100,000 y
        
        """
        
        wastemass = 0.0
        activity = {1e2: 0.0, 1e5: 0.0}
        strm = self.discharge_streams()
        strm = self.cool(strm)
        if(self.data['reprocessing']!='none'):
            waste = self.reprocess(strm, **kwargs)
        else:
            waste = strm
        # sum all streams in waste to get one big stream; this decreases
        # computation time relative to running separate ORIGEN cases for each
        waste = self._sum_stage_waste_streams(waste)
        # waste radioactivity
        for time in activity: # all fuel cycles should be cooled 5 y
            decaytime = time-float(waste.comp.columns[-1])
            temp = self.decay(waste,endtime=decaytime,steps=np.log10(time)*10)
            timestring = min(temp.act.columns, \
                             key=lambda x:abs(float(x)-time))
            activity[time] = sum(temp.act[timestring])
        return(waste.mass, Counter(activity))
    
    
    def _sum_stage_waste_streams(self, str_list, time=5.0):
        """Sum the streams in a list (or list of lists)"""
        
        # time is preset because all fces eg have cooling time of 5y
        stg = self.evaluationgroup+'-'+str(self.number)
        strm = stream.empty()
        if(not hasattr(str_list, 'index')): # # recast as list; quacks?
            str_list = [str_list]
        for item in str_list:
            if(isinstance(item, Sequence)):
                strm = strm.add(self._sum_stage_waste_streams(item))
            else:
                strm = strm.add(item)
        return(strm)


    #--------------------------------------------------------------------------
    # DISCHARGE
    # Given the hidden stream data contained in the stage instance hidden data
    # dictionary, this function returns the isotopic composition for any
    # material streams output from the irradiation system in that stage.
    #--------------------------------------------------------------------------
    def discharge_streams(self):
        """Get data for each stream discharged from the irradiation system
        
        Parameters
        ----------
        (None)
        
        Returns
        -------
        Stream instance (or list of Stream instances) discharged 
            from stage irradiation system
        
        """
        
        return(self.reactor.discharge_streams())
    
    
    #--------------------------------------------------------------------------
    # ORIGEN INTERFACE
    # This Python package uses ORIGEN to calculate time-dependent properties
    # of material streams. This is important for two applications:
    #
    # 1. Cooling between reactor discharge and separations/waste loading
    # 2. Radioactive decay in the waste after waste form preparation
    #
    # Both of these are accomplished using the decay method below. Application
    # (2) uses "decay" directly, while application (1) uses a special method,
    # "cool", which uses specific inputs from the stage data. These methods
    # make the ORIGEN input file, run ORIGEN, and extract the data to produce
    # new streams.
    #--------------------------------------------------------------------------
    def decay(self, strm, starttime='last', endtime=500.0,
              steps=10, interp='L', **kwargs):
        """Calculate waste properties over periods of time for streams;
        streams can be held in lists of any arrangement and depth, and this
        function will return them in the same organization.
        
        Parameters
        ----------
        strm: list (or Sequence type) OR single Stream instance
        
        starttime (optional): float or str
            Represents column header in pandas dataframe attribute from
            which to pull composition data for ORIGEN calculation.
        
        endtime (optional): float
            Final time for the decay calculation
    
        steps (optional): float/int
            Number of steps required for the decay calculation; subject to
            ORIGEN requirements as well
            
        interp (optional): str
            Indicates interpolation method between decay calculation steps
        
        **kwargs
        --------
        rerun (optional): Boolean
            Indicates whether ORIGEN should be called or whether the
            stream can be updated with existing results
        
        cooling (optional): Boolean
            Flag to pass to origen.py module which controls the generation
            of the ORIGEN input file
        
        Returns
        -------
        A single Stream instance, or a list of Stream instances, with data:
        - composition
        - decay heat
        - radioactivity
        - radiotoxicity
        
        """
        
        if(not hasattr(strm, 'index')): # recast as list; quacks?
            strm = [strm]
        out = []
        for item in strm:
            if(isinstance(item, Sequence)):
                out.append(self.decay(item, starttime, endtime,
                                      steps, interp, **kwargs))
            else:
                out.append(self._decay_stream(item, starttime, endtime,
                                              steps, interp, **kwargs))
        if(len(out)==1): # return the output streams
            return(out[0])
        else:
            return(out)


    def cool(self, str_list, starttime='last', rerun=True):
        """Run cooling cases for streams; streams can be held in lists of any
        arrangement and depth, and this function will return them in the same
        organization.
        
        Parameters
        ----------
        str_list: list (or Sequence type) OR single Stream instance
            One or more Streams to be cooled
        
        starttime (optional): float or str
            Represents column header in pandas dataframe attribute from
            which to pull composition data for ORIGEN calculation.
        
        rerun (optional): Boolean
            Indicates whether ORIGEN should be called or whether the
            stream can be updated with existing results

        Returns
        -------
        A single Stream instance, or a list of Stream instances
        
        """
            
        ct = self.data['coolingtime']
        return(self.decay(str_list, starttime=starttime, endtime=ct,
                          steps=ct, interp='I', rerun=rerun, cooling=True,
                          remove=True))


    def _decay_stream(self, strm, start, end, steps, interp, **kwargs):
        """Run decay calculation for single stream"""
        
        assert end >= 0.0, 'Time must be greater than or equal to zero.'
        assert steps < 999, 'ORIGEN allows a maximum of 999 time steps.'
        if(start == 'last'):
            start = strm.comp.columns[-1]
        if(kwargs.get('cooling')):
            cool=True
        else:
            cool=False
        inp = origen.Origen.make_file(self.wdpath, self.name, strm, start,
                                      end, steps, interp, cooling=cool)
        if(not kwargs.get('rerun')==False):
            if(kwargs.get('remove')):
                self._remove_old_outfiles(inp)
            origen.Origen.run_origen(self.wdpath, inp)
        new_strm = origen.Origen.update_stream(self.wdpath, strm, inp, start)
        return(new_strm)

    
    def _remove_old_outfiles(self, infile):
        """Ensure no ORIGEN in/out files match the new case to be run"""
        
        file = infile[:-4]
        files_in_outdir = [x for x in os.listdir(self.wdpath) if file in x]
        files_in_outdir = [x for x in files_in_outdir if '.inp' not in x]
        if(len(files_in_outdir) > 0):
            for i in range(0, len(files_in_outdir)):
                #print files_in_outdir[i]
                os.remove(os.path.join(self.wdpath, files_in_outdir[i]))


    #--------------------------------------------------------------------------
    # REPROCESSING
    # In some of the stages, material streams output from the irradiation
    # system undergo separations in which some of the materials are recovered,
    # altering the mass and composition of the waste. The following functions
    # control the application of separations to material streams based on
    # separation data. The data is in the form of fractions of the material
    # stream that end up in each waste stream. The data is held in text files
    # in the data/sep/ directory, named by the type of separation and the
    # recovered species.
    #--------------------------------------------------------------------------
    def reprocess(self, str_list, **kwargs):
        """Create dict(s) of isotopic composition, heat for unrecovered 
        materials in stream(s) that result from separations applied in the 
        stage.
        
        Parameters
        ----------
        str_list: Stream instance, or list of Stream instances
            List of streams required for the MSR stages; 
            single Stream instance acceptable for other cases.
        
        kwargs:
        include: str
            Indicate whether the reprocess method should return just the waste
            streams (default) or, for some separations processes, just the 
            products (="products") or both wastes and products (="all").
        
        Returns
        -------
        A single Stream instance, or a list of Stream instances
        
        """
        
        def flatten(x): # flatten the lists returned from this method
            if isinstance(x, Sequence):
                return [a for i in x for a in flatten(i)]
            else:
                return [x]
        ostrms = [] # output streams
        #if(self.data['reactor']=='msr'): # molten salt special case
        #    try: # requires [df1, fp2]
        #        ostrms.append(self.separation.reprocess(str_list))
        #    except: # list passed as [fp1, [df1, fp2]]
        #        ostrms.append(str_list[0]) # discharged
        #        ostrms.append(self.separation.reprocess(str_list[1]))
        #else:
        if(not hasattr(str_list, 'index')): # recast as list
            str_list = [str_list]
        for item in str_list:
            if(isinstance(item, Sequence)):
                ostrms.append(self.reprocess(item, **kwargs))
            else:
                if(not any (x in item.form for x in ['df', 'fp'])):
                    # df streams  not reprocessed
                    # expected catches: DF from BNBSFR, MSR;
                    # MSR salt treatment FP
                    try:
                        ostrms.append(self.separation.reprocess(item,**kwargs))
                    except AttributeError: # no self.separation
                        ostrms.append(item) # do not reprocess stream
                else: # df goes to ostrms unmodified
                    ostrms.append(item)
        if(len(ostrms)==1):
            return(ostrms[0])
        else:
            return(flatten(ostrms))

    
    #--------------------------------------------------------------------------
    # WASTE LOADING
    # After discharge and separations, waste must be loaded into canisters for
    # disposal. The amount of waste loaded per canister governs the properties
    # of the canister after disposal (heat generation rate, radionuclide
    # source terms, etc.). Different spent fuels will be loaded in different
    # amounts; different reprocessing methods will produce waste streams with
    # varying characteristics that influence the choice of waste form and the
    # amount of waste nuclides that it can contain. Based on these factors,
    # the amount of waste per canister is calculated.
    #--------------------------------------------------------------------------
    def load_waste(self, str_list, time='last', **kwargs):
        """Load a list of waste streams (SNF/HLW), one by one, into a 
        wasteform suitable for long-term disposal in a geological repository.
        
        Note that if multiple streams are input, any keyword args will be
        applied to all of them. Therefore, it's recommended that this method
        be called on single streams rather than lists.

        Paramters
        ---------
        str_list: list or tuple
            List of Stream instances
        
        time (optional): str or float
            Time to grab Stream instance data to transfer to waste form
            (must correspond to dataframe column or keyword)
        
        kwargs
        --------
        verbose: Boolean
            Print statements during calculation
        
        loading (for SNF): float or int or str
            If number: number of assemblies per canister; must agree with
            value in data file unless kwarg 'can' is also specified.
        
        can (for SNF): dict
            Keys 'Diameter' and 'Length'; both values floats
        
        consolidate (for HTGR SNF): Boolean
            If True and reactor is HTGR, consolidate fuel particles from
            prismatic blocks to increase loading.
        
        plot (for aqueous glass and e-chem metal HLW): Boolean
            If True, generates plot of loading solution space

        recycle (for e-chem ceramic HLW): Boolean
            If True, indicates distillation and reclamation of carrier salt, 
            and therefore concentration of the waste.
            
        Returns
        -------
        List of waste form instance with attributes on a per-canister basis
        
        """
        
        if(not hasattr(str_list, 'index')): # recast as list; quacks?
            str_list = [str_list]
        if(self.data['reactor'] in ['msr','ffh']
           and 'salttreatment' in self.data.keys()):
            outstreams = self._load_msr_wastes(str_list, time, **kwargs)
        else:
            outstreams = []
            for item in str_list:
                if(isinstance(item, Sequence)):
                    outstreams.append(self.load_waste(item, **kwargs))
                else:
                    outstreams.append(self._load_stream(item, time, **kwargs))
        # remove None elements from outstreams
        outstreams = [x for x in outstreams if x is not None]
        if(len(outstreams)==1): # return the output streams
            return(outstreams[0])        
        else:
            return(outstreams)


    def _load_stream(self, str_inst, time, **kwargs):
        """Load waste stream (SNF/HLW) into a wasteform suitable for long-term
        disposal in a geological repository.
        
        Paramters
        ---------
        str_inst: Stream instance
        
        time: str or float
            Time to grab Stream instance data to transfer to waste form
            (must correspond to dataframe column or keyword)
        
        **kwargs:
        verbose (optional): Boolean
            Print statements during calculation
        
        plot (optional): Boolean
            Plot the loading constraints and solution region
        
        loading_level (optional): str
            'high' or 'low' loading for HTGR SNF; default 'low'
        
        recycle (optional): Boolean
            For ceramic wastes, 'True' indicates distillation and reclamation 
            of carrier salt, and therefore concentration of the waste.
        
        Returns
        -------
        Waste form instance with attributes on a per-canister basis
        
        """

        # SNF CASE
        if(str_inst.id=='snf' and (self.data['reprocessing']=='none' or
                                   self.data['reactor']=='sfr')):
            wl = loading.Loading(data=self.data, datapath=self.datapath)
            if(self.data['reactor']=='htgr' and
               'loading_level' not in kwargs.keys()):
               kwargs['loading_level'] = 'low'
        # HLW CASES - waste forms are differentiated by reprocessing method
        else:
            # UREX and other aqueous methods
            if(self.data['reprocessing'] in ['urex','urex+','thorex']
               and 'glass' in str_inst.form):
                wl = loading.AqGlass(data=self.data, datapath=self.datapath)
            # Electrochemical processing
            elif(self.data['reprocessing']=='echem'
                 and 'metal' in str_inst.form):
                wl = loading.EcMetal(data=self.data, datapath=self.datapath)
            elif(self.data['reprocessing']=='echem'
                 and 'ceramic' in str_inst.form):
                wl = loading.EcCeramic(data=self.data,datapath=self.datapath)
            # Melt-refining
            elif(self.data['reprocessing']=='meltrefining'
                 and 'skull' in str_inst.form):
                wl = loading.Skull(data=self.data,datapath=self.datapath)
            elif(self.data['reprocessing']=='meltrefining'
                 and 'gas' in str_inst.form):
                wl = loading.CapturedCs(data=self.data,datapath=self.datapath)
            # Unsupported methods
            else:
                #raise NotImplementedError('Loading for wastes not supported')
                print('Loading for '+str_inst.form+' not supported')
                return # nothing
        return(wl.load_waste(str_inst, time, **kwargs))
        

    def _load_msr_wastes(self, str_list, time, **kwargs):
        """Load MSR wastes
            
        Parameters
        ----------
        str_list: flat list of all MSR waste streams to be loaded
        
        time: str or float
            Time to grab Stream instance data to transfer to waste form
            (must correspond to dataframe column or keyword)
        
        **kwargs:
        - loading_fraction: float
            Specify the HLW loading fraction
        
        Returns
        -------
        Waste form instances with attributes on a per-canister basis

        """
        
        # salt treatment wastes
        fp1 = [x for x in str_list if x.form=='fp1'][0]
        salttreat = loading.MSRMetal(data=self.data, datapath=self.datapath)
        fp1_wf = salttreat.load_waste(fp1, **kwargs)
        # ceramic wastes - distill fuel salt?
        c_list = [x for x in str_list if x.form!='fp1']
        ceramic = c_list[0]
        if(len(c_list) > 1): # should have max length of two
            ceramic = ceramic.add(c_list[1], time=time)
        loadceramic=loading.MSRCeramic(data=self.data,datapath=self.datapath)
        ceramic_wf = loadceramic.load_waste(ceramic, **kwargs)
        return(fp1_wf, ceramic_wf)


###############################################################################
