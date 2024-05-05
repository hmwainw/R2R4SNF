###############################################################################
"""
    Last edited on September 15, 2020

    @author: matz

    comments: Waste loading classes and methods
    
"""
###############################################################################
import os
import imp
import numpy as np
import pandas as pd
import scipy
from scipy import optimize
from copy import deepcopy
from nwpy import nuclidedata
from nwpy import plot
from nwpy import stream
###############################################################################


class Loading(object):
    
    
    """
    Before long-term storage and final disposal, waste must be loaded into
    canisters. For SNF, this means some number of fuel assemblies will be
    put into each canister. For HLW from reprocessing, this means that some
    amount of waste, mixed with a matrix to produce a suitable waste form, is
    poured into each canister.
    
    This class contains the methods and data required to calculate waste 
    loading. These methods are unique with respect to the condition of the
    stream to be loaded. For example, the calculation of waste loading for
    the aqueous HLW from UREX reprocessing is different than that for the
    ceramic and metal HLW from electrochemical reprocessing. However, all 
    are controlled by a central method called loadWaste. The Base class 
    considers the loading of spent fuel, but special classes for other
    wastes can be inherited if necessary.
    
    Given knowledge about the waste stream and desired waste form, this
    general class inherits methods associated with one of the specific
    waste loading classes below. Because knowledge of the waste stream is
    required, this class is not instantiated in the Stage instance until
    it is required there.

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    NOTE: The number of canisters will almost *never* be an integer. The
    number of canisters is determined by first calculating how much material
    can go into a single canister, then dividing the total mass to be loaded
    by that value.
    
    Because fractions of canisters are not realistic, the total canister
    count should be rounded up to the nearest integer. However, this 
    presents an issue, as one of two things must happen:
    (1) if mass is conserved and canisters are assumed to all be the same, 
        the amount of material in each canister must be *lower* than the 
        value previously calculated;
    (2) if mass is conserved but the canisters are not all the same, the 
        extra material can be placed into one canister loaded to a lesser
        extent than all the others
    (3) Mass is not conserved but all canisters are the same and because 
        there is fractionally more canisters than calculated, the total 
        mass of material distributed across those canisters increases.
        
    Of these options, option (1) is the best option. Option (2) presents
    difficulties in later calculations, and option (3) may significantly
    bias the results, especially if only a few canisters are required 
    such that the fractional increase of material mass is significant. 
    The result of option (1) is the violation of the assumptions and 
    input data that went into performing the calculation. However, this 
    method offers more accurate results when the canister count is small, 
    and when the canister count is large the effect should be minimal 
    because the fractional increase in canister count is insignificant 
    relative to the total number of canisters. Additionally, the input
    data and assumptions are themselves subject to distributions and 
    uncertainty - small changes in their values will not be unrealistic 
    and should not affect the meaningfulness of the conclusions.

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    """
    
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
    

    def load_waste(self, str_inst, time='last', **kwargs):
        """Load SNF into canisters
        
        Parameters
        ----------
        str_inst: Stream instance
            SNF stream to be loaded into canister
        
        time (optional): str or float
            Time to grab Stream instance data to transfer to waste form
            (must correspond to dataframe column or keyword)
        
        **kwargs
        --------
        loading: float or int or str
            If number: number of assemblies per canister; must agree with
            value in data file unless kwarg 'can' is also specified.
            
        can: dict
            Keys 'Diameter' and 'Length'; both values floats
        
        consolidate: Boolean
            If True and reactor is HTGR, consolidate fuel particles from
            prismatic blocks to increase loading.
            
        verbose (optional): Boolean
            Print statements during calculation
        
        Returns
        -------
        WasteForm instance
        
        """
        
        loaddata = self._get_loading_data(str_inst.id, str_inst.form)
        if(not kwargs.get('loading')): # use default value if not specied
            kwargs['loading'] = loaddata.asm_per_canister
        if(kwargs.get('can')): # get the can specs based on input/data
            self.can = kwargs['can']
        else:
            self.can = self.can[int(kwargs['loading'])]
        if(self.data['reactor']=='htgr' and kwargs.get('consolidate ')):
            v_can = np.pi*((self.can['Diameter']/2)-self.can['Thickness'])**2
            v_can *= self.can['Length']-self.can['Thickness']*2
            hm_loading = loaddata.fuel_dens*v_can
        else:
            hm_loading = kwargs['loading']*loaddata.hm_per_asm
        n_canisters = str_inst.mass/hm_loading
        return(self._make_wasteform(str_inst, n_canisters, time, **kwargs))
    
        
    def _get_loading_data(self, streamid, streamform):
        """Get the loading data based on the stream and stage data"""
        
        if(any([x in streamform for x in ['driver', 'blanket']])):
            streamform = self._remove_dfbb_modifier(streamform)
        dp = os.path.join(self.datapath,'load')
        if(self.data['reactor'] in ['msr','ffh']
           and 'salttreatment' in self.data.keys()):
            if(streamform=='fp1'):
                #file = self.data['reprocessing']+'_metal.py'
                file = 'msr_metal.py'
            else:
                #file = self.data['reprocessing']+'_ceramic.py'
                file = 'msr_ceramic.py'
        elif(streamid=='snf'):
            file = 'snf_'+self.data['reactor']+'.py'
        else:
            file = self.data['reprocessing']+'_'+streamform+'.py'
        data = imp.load_source(file, os.path.join(dp,file))
        try:
            self.can = data.canister
        except:
            self.can = None
        return(data)
    
    
    @staticmethod
    def _remove_dfbb_modifier(streamform):
        """Identify form of stream by stripping driver/blanket modifier"""
    
        dfbb_id = [s for s in ['blanket', 'driver'] if s in streamform][0]
        tag = streamform.find(dfbb_id)
        return(streamform[:tag-1]) # also remove the underscore
    

    def _make_wasteform(self, strm, n_can, time, **kwargs):
        """Produce a WasteForm instance by distributing the waste Stream 
        data over the calculated number of canisters"""
        
        n_can = np.ceil(n_can) # make sure n_can is int
        mass = strm.mass/n_can # calculate mass per canister
        if(kwargs.get('verbose')):
            self._print_result(mass, n_can)
        wf = stream.WasteForm(mass,n_can,self.can,id=strm.id,form=strm.form,
                              evaluationgroup=strm.evaluationgroup,
                              stagenumber=strm.stagenumber)
        if(kwargs.get('loading')): # number of wasteforms (prob SNF) per can
            wf.loading = kwargs['loading']
        if(kwargs.get('loading_fraction')): # fractional loading of HLW
            wf.loading_fraction = kwargs['loading_fraction']
        if(hasattr(strm, 'batch')): # batch number for sfr bnb wfs
            wf.batch = strm.batch
        series = stream.get_srs(strm, time) # modify, reassign data to wf
        for attr, old_srs in series:
            new_srs = old_srs/n_can
            setattr(wf, attr, new_srs.to_frame())
        return(wf)
    
    
    @staticmethod
    def _print_result(m, n):
        """Print the results of the waste loading optimization calculation;
        
        Parameters
        ----------
        m: mass per canister
        n: number of canisters
            
        """
        
        print('Max waste per canister (kg): '+str(round(m/1e3,4)))
        print('HLW Canisters: '+ str(int(np.ceil(n))))
    

    def easy_hlw_loading(self, str_inst, lf):
        """Given user-specified loading fraction, calculate the number of
        HLW waste packages - for HLW ONLY
        
        Parameters
        ----------
        str_inst: Stream instance
        
        lf: float
            Loading fraction
            
        """
        
        assert lf > 0.0, 'User-supplied loading fraction must be greater than 0'
        assert lf < 1.0, 'User-supplied loading fraction must be less than 1'
        loaddata = self._get_loading_data(str_inst.id, str_inst.form)
        mass_wf = str_inst.mass/lf
        n_can = mass_wf/loaddata.canister['Mass limit']
        return(n_can)


################################################################################
# GLASS from UREX separation
################################################################################


class AqGlass(Loading):


    """
    Glass HLW from aqueous reprocessing is loaded with fission product oxides
    according to constraints in the load datafile.
    
    """


    def __init__(self, **kwargs):
        super(AqGlass, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)
    
    
    # Glass from aqueous reprocessing is loaded w/ FP oxides according
    # to constraints in the loading.py file in the data directory.
    def load_waste(self, str_inst, time='last', **kwargs):
        """Formulate the constraint functions that bound
        glass waste loading (from urex)
        
        Parameters
        ----------
        self: AqGlass loading instance
        
        str_inst: Stream instance
            Waste from aqueous separations to be loaded into glass form
            
        time (optional): str or float
            Time to grab Stream instance data to transfer to waste form
            (must correspond to dataframe column or keyword)
        
        **kwargs (optional):
        - verbose: if True, print details about waste loading
        - plot: if True, generates plot of loading solution space
        - loading_fraction: float
            User-specified waste loading fraction
            
        """

        if(kwargs.get('loading_fraction')): # ignore linear optimization
            n_can = self.easy_hlw_loading(str_inst, kwargs['loading_fraction'])
        else:
            # Isolate important composition, heat series
            comp_srs = str_inst.comp[stream._get_time(time, str_inst.comp)]
            heat_srs = str_inst.heat[stream._get_time(time, str_inst.heat)]
            # Import required data
            A, b, key, oxide_stoich = self._import_constraints(str_inst)
            # Make oxides
            ox_srs = self._oxidize(comp_srs, oxide_stoich)
            # Update the constraint matrix structures
            A, b = self._update_constraints(ox_srs, heat_srs, A, b)
            res, n_can = self._evaluate_constraints(A, b, ox_srs)
            if(kwargs.get('verbose')): # additional printouts for this method
                print('Max possible loading: '+str(round(res.x[0]/sum(res.x),4)))
                print('Max oxide per canister (kg): '+str(round(res.x[0]/1e3,4)))
            if('plot' in kwargs.keys()): # plot if asked
                plot.LoadPlot.plot(A, b, key, units='grams')
            #matrix_mass = res.x[1] #matrix_mass=2900.0-res.x[0]
        return(self._make_wasteform(str_inst, n_can, time, **kwargs))
    

    def _import_constraints(self, str_inst):
        """Import systems of equations, and necessary loading data"""
        
        loaddata = self._get_loading_data(str_inst.id, str_inst.form)
        A_ub = loaddata.glass['ineq'][:,:2]
        b_ub = loaddata.glass['ineq'][:,2]
        key_ub = loaddata.glass['key']
        return(A_ub, b_ub, key_ub, loaddata.oxide)
    
    
    def _update_constraints(self, oxide_srs, heat_srs, A_ub, b_ub):
        """Update system of equations w stream-specific information"""

        # Entries 1-5 in the constraints already done; only 6-8 need work
        # Heat (totalheat/oxidemass)
        A_ub[4][0]=sum(heat_srs)/sum(oxide_srs)
        # MoO3 fraction
        try:
            x_moo3 = oxide_srs['moo3']/sum(oxide_srs)
        except: # no moly in glass
            x_moo3 = 0.0
        A_ub[5] = [(x_moo3-b_ub[5]), -b_ub[5]]
        b_ub[5] = 0.0
        # Rare earth fraction
        try:
            m_reox = (oxide_srs['rh2o3']+oxide_srs['ruo2']+
                      oxide_srs['ag2o']+oxide_srs['pdo'])
        except: # no rare earth in glass
            m_reox = 0.0
        x_reox = m_reox/sum(oxide_srs)
        A_ub[6] = [(x_reox-b_ub[6]), -b_ub[6]]
        b_ub[6] = 0.0
        return(A_ub, b_ub)
    
    
    @staticmethod
    def _evaluate_constraints(A_ub, b_ub, oxide_srs):
        """Solve system of linear equations to get the number of canister
        required for the stream and the mass loaded in each canister."""
        
        res = scipy.optimize.linprog([-1, 0], A_ub=A_ub, b_ub=b_ub)
        n_canisters = sum(oxide_srs)/res.x[0]
        return(res, n_canisters)
    
    
    # The radionuclides separated in aqueous reprocessing are loaded into the
    # glass waste form as oxides; because the oxygen adds diluent mass to the
    # stream, the new mass of the oxidized stream must be calculated.
    def _oxidize(self, comp, stoich):
        """Calculate the masses of waste stream oxides
        
        Parameters
        ----------
        self: Loading instance
        
        comp: Pandas Series
            Composition of the Stream to be oxidized
            
        stoich: dict
            Stoichiometry of elemental oxides
        
        """

        ox = {}
        for nuc in comp.index:
            el, a, meta = nuclidedata.determine_nuclide_info(nuc)
            if(comp[nuc]==0): # avoid gases, etc
                continue
            else:
                multiplier = ((stoich[el][0]*float(a)+stoich[el][1]*16)/
                              float(a)/stoich[el][0])
            key = self._build_oxide(el, stoich[el])
            try: # add to dictionary
                ox[key]+=comp[nuc]*multiplier
            except KeyError:
                ox[key]=comp[nuc]*multiplier
        srs = pd.Series(list(ox.values()),index=list(ox.keys()),name=comp.name)
        return(srs)
    
    
    @staticmethod
    def _build_oxide(element, stoich):
        """Build oxide chemical symbol"""
        
        ox = element
        if(stoich[0]!=1):
            ox += str(stoich[0])
        ox += 'o'
        if(stoich[1]!=1):
            ox += str(stoich[1])
        return(ox)


################################################################################
# METAL from electrochemical reprocessing
################################################################################


class EcMetal(Loading):
    
    
    """
    Metal HLW from electrochemical reprocessing is loaded with noble metal 
    fission products, along with Zr used in the alloying of metal fuels and
    HT9 cladding, according to constraints in the load datafile.
    
    """
    
    
    def __init__(self, **kwargs):
        super(EcMetal, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)


    def load_waste(self, str_inst, time='last', **kwargs):
        """Formulate the constraint functions that bound
        metal waste loading (from pyroprocess)
        
        Parameters
        ----------
        self: EcMetal loading instance
        
        str_inst: Stream instance
            Electrochemical separations waste to be loaded into metal form
        
        time (optional): str or float
            Time to grab Stream instance data to transfer to waste form
            (must correspond to dataframe column or keyword)
        
        **kwargs (optional):
        - verbose: if True, print details about waste loading
        - plot: if True, generates plot of loading solution space
        - loading_fraction: float
            Use user-specified loading fraction
        
        """
        
        if(kwargs.get('loading_fraction')): # ignore cladding, hardware
            n_can = self.easy_hlw_loading(str_inst, kwargs['loading_fraction'])
        else:
            # Isolate important composition series
            comp_srs = str_inst.comp[stream._get_time(time, str_inst.comp)]
            # Get linear data structures, update based on stream data
            Au, bu, keyu, Ae, be, keye, m, x = self._import_constraints(str_inst)
            Ae = self._update_constraints(comp_srs, Ae, m, x)
            res, n_can = self._evaluate_constraints(Au, bu, Ae, be, m)
            if(kwargs.get('plot')): # plot if asked
                plot.LoadPlot.plot(Au, bu, keyu, Ae, be, keye)
        return(self._make_wasteform(str_inst, n_can, time, **kwargs))

    
    def _import_constraints(self, str_inst):
        """Import systems of equations, and necessary loading data"""
        
        loaddata = self._get_loading_data(str_inst.id, str_inst.form)
        A_ub = loaddata.metal['ineq'][:,:2]
        b_ub = loaddata.metal['ineq'][:,2]
        key_ub = loaddata.metal['key']
        A_eq = loaddata.metal['eq'][:,:2]
        b_eq = loaddata.metal['eq'][:,2]
        key_eq = loaddata.metal['key'][len(A_ub):]
        # masses of key components;
        # Heavy metal in Stage (t -> kg) is self.data['mass']*1e3
        m = {'nmfp': str_inst.mass, # g Noble metal FP in stream
             'zr': self.data['mass']*loaddata.zr_to_hm, # Alloy Zr; not FP
             'ht9': self.data['mass']*loaddata.ht9_to_hm # cladding
             }
        x = deepcopy(m)
        x.update((k,v/sum(m.values())) for k,v in x.items())
        return(A_ub, b_ub, key_ub, A_eq, b_eq, key_eq, m, x)
    
    
    @staticmethod
    def _update_constraints(comp, A_eq, m, x):
        """Update system of equations w stream-specific information"""
        
        # fraction of zr that is nmfp (mass in kg)
        mass_zrfp = sum([comp[k] for k in comp.index if 'zr' in k])
        x_zrfp = mass_zrfp/sum(m.values())
        A_eq[0] = [1, -0.15/(x['zr']+x_zrfp-0.15)]
        return(A_eq)
    
    
    @staticmethod
    def _evaluate_constraints(A_ub, b_ub, A_eq, b_eq, m):
        """Solve system of linear equations"""
        
        # solve the system
        res = scipy.optimize.linprog([-1, 0], A_ub=A_ub, b_ub=b_ub,
                                     A_eq=A_eq, b_eq=b_eq) # kg
        try:
            n_canisters = int(np.ceil(sum(m.values())/res.x[0])) # kg
        except OverflowError:
            print(res)
            print(m)
            raise
        return(res, n_canisters)


################################################################################
# CERAMIC from electrochemical reprocessing
################################################################################


class EcCeramic(Loading):
    
    
    """
    Ceramic HLW from electrochemical reprocessing is loaded with active
    metal fission products that dissolve into the electrorefiner salt,
    along with some actinides that are not recovered.
    
    """
    
    
    def __init__(self, **kwargs):
        super(EcCeramic, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)


    def load_waste(self, str_inst, time='last', **kwargs):
        """Calculate waste loading in sodalite ceramic (from pyroprocess)
        
        Parameters
        ----------
        self: AqGlass loading instance
            
        str_inst: Stream instance
            Electrochemical separations waste to be loaded into ceramic
        
        time (optional): str or float
            Time to grab Stream instance data to transfer to waste form
            (must correspond to dataframe column or keyword)
            
        **kwargs (optional):
        - verbose: if True, print details about waste loading
        - recycle: if True, recycle salt, increasing loading level
        - loading_fraction: float
            User-specified waste loading fraction
        """

        loaddata = self._get_loading_data(str_inst.id, str_inst.form)
        if(kwargs.get('loading_fraction')):
            n_can = self.easy_hlw_loading(str_inst, kwargs['loading_fraction'])
        else:
            #if(kwargs.get('recycle')):
            #    x_fp = loaddata.sodalite['salt high fp fraction']
            #else:
            #    x_fp = loaddata.sodalite['salt fp fraction']
            x_fp = loaddata.sodalite['salt high fp fraction']
            x_slz = loaddata.sodalite['SLZ salt fraction']
            x_cwf = loaddata.sodalite['CWF zeolite fraction']
            kwargs['loading_fraction']=x_fp*x_slz*x_cwf
            mass_cwf = str_inst.mass/kwargs['loading_fraction']
            n_can = mass_cwf/loaddata.sodalite['Canister mass limit']
        return(self._make_wasteform(str_inst, n_can, time, **kwargs))


################################################################################
# CERAMIC from Molten Salt Reactor online reprocessing
################################################################################


class MSRCeramic(Loading):
    
    
    """
    Ceramic HLW in fluorapatite made from the FP recovered from MSR salt. 
    
    In the MSR separations process, a portion of the circulating salt is
    diverted from the primary loop to undergo fission product removal. 
    Uranium is recovered from the salt by fluorination, after which the 
    barren salt undergoes distillation in which rare earth FP are left at 
    the still bottoms. The distilled salt is recombined with uranium and
    reintroduced to the reactor; the FP are left to accumulate in waste
    salt in a tank adjacent to the column.
    
    Any fuel salt discarded as spent fuel is assumed to join the fission 
    product-loaded salt in the accumulator.
    
    Before loading, the accumulator salt will likely be recovered (perhaps
    by distillation), leaving only FP and actinides. These waste species 
    will be embedded in some robust waste form - in this case that waste 
    form is assumed to be a fluorapatite material, due to expected chemical 
    compatability with fission product fluorides and any residual salt.
        
    """
    
    
    def __init__(self, **kwargs):
        super(MSRCeramic, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)


    def load_waste(self, str_inst, time='last', **kwargs):
        """Calculate waste loading in fluorapatite ceramic (from MSR)
        
        Parameters
        ----------
        str_inst: Stream instance
            Fission products from MSR separations AND/OR discharged fuel salt
        
        time (optional): str or float
            Time to grab Stream instance data to transfer to waste form
            (must correspond to dataframe column or keyword)
        
        **kwargs (optional):
        - verbose: if True, print details about waste loading
        - loading_fraction: float
            User-specified waste loading fraction
        """
        
        temp = deepcopy(str_inst)
        temp.id = 'hlw'
        if(kwargs.get('loading_fraction')):
            n_can = self.easy_hlw_loading(temp, kwargs['loading_fraction'])
        else:
            temp = deepcopy(str_inst)
            loaddata = self._get_loading_data(temp.id, temp.form)
            x_w = 1.0 # no excess salt in final waste stream
            m_w = temp.mass/x_w # kg
            m_fluorapatite = m_w/loaddata.fluorapatite['FP loading']
            n_can = m_fluorapatite/loaddata.fluorapatite['Canister mass limit']
        if('verbose' in kwargs.keys()):
            #fp_load = temp.mass_fraction('fp')*temp.mass/m_fluorapatite
            print('MSR separations and discharged salt waste loading')
            #print('Maximum FP loading: '+str(round(fp_load,4)))
            print('HLW canisters: '+str(int(np.ceil(n_can))))
        # scale masses, heats to a per canister basis
        temp.form = 'ceramic'
        return(self._make_wasteform(temp, n_can, time, **kwargs))


################################################################################
# WASTE FROM MSR SALT TREATMENT
################################################################################


class MSRMetal(Loading):
    
    
    """
    As the fuel salt circulates the MSR core, fission products that do not
    form stable fluorides in the salt must be removed. Within the primary
    loop, salt is treated to continuously remove fission products that are
    not stable in the salt. Noble gases are sparged from the fuel salt by
    bubbling helium and are held in a tank; noble metal fission products 
    plate out of the salt and are disposed of in a metal waste form, which
    is described here.
    
    """
    
    
    def __init__(self, **kwargs):
        super(MSRMetal, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)


    def load_waste(self, str_inst, time='last', **kwargs):
        """Calculate waste loading in metal waste form; gases are sent to
        decay tanks (not considered, metals go to waste form
        
        Parameters
        ----------
        str_inst: Stream instance
            Fission products separated from online salt treatment
        
        time (optional): str or float
            Time to grab Stream instance data to transfer to waste form
            (must correspond to dataframe column or keyword)
        
        **kwargs (optional):
        - verbose: if True, print details about waste loading
        - loading_fraction: float
            User-specified waste loading fraction
        """
        
        temp = deepcopy(str_inst)
        series = stream.get_srs(temp, time)
        # separate gases
        gases = ['h', 'he', 'ne', 'ar', 'kr', 'xe', 'rn']
        to_drop = []
        for attr, srs in series:
            for nuc in srs.index:
                el = nuclidedata.determine_nuclide_info(nuc)[0]
                if(el in gases):
                    to_drop.append(nuc)
            srs = srs.drop(to_drop) # to the tank
            setattr(temp, attr, srs.to_frame())
        temp.mass = sum(temp.comp[stream._get_time(time, temp.comp)])
        if(kwargs.get('loading_fraction')):
            n_can = self.easy_hlw_loading(temp, kwargs['loading_fraction'])
        else:
            loaddata = self._get_loading_data(temp.id, temp.form)
            m_metal = temp.mass/loaddata.metal['FP loading']
            n_can = m_metal/loaddata.metal['Canister mass limit']
        if(kwargs.get('verbose')):
            #fp_load = temp_str.mass_fraction('fp')*temp.mass/m_metal
            print('MSR Salt Treatment Metal FP Loading')
            #print('FP loading: '+str(round(fp_load,4)))
            print('HLW canisters: '+str(int(np.ceil(n_can))))
        # scale masses, heats to a per canister basis
        temp.form = 'metal'
        return(self._make_wasteform(temp, n_can, time, **kwargs))


################################################################################
# CAPTURED CS+RB OFF-GAS WASTE FROM MELT REFINING
################################################################################


class CapturedCs(Loading):


    """
    In melt refining, the volatile and semi-volatile FP are released from the
    melted fuel as gases. There are four main groups of volatile FP: alkali
    metals (Rb, Cs), halogens (I, Br), noble gases (Kr, Xe), and cadmium.
    The noble gases are held in tanks to allow for decay and, ultimately,
    controlled release to the environment. The remaining radionuclides will 
    be stripped from the effluent. Of these, the alkali elements have by
    far the largest activity and decay heat. Therefore, the consideration 
    of waste forms from the melt off gas consider the capture of Cs (and Rb) 
    onto molecular sieves and, assuming these sieves to be an acceptable waste
    form, the subsequent emplacement of those loaded sieves into canisters 
    for disposal.
    
    """
        
        
    def __init__(self, **kwargs):
        super(CapturedCs, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)


    def load_waste(self, str_inst, time='last', **kwargs):
        """Calculate waste loading in molecular sieve
        
        Parameters
        ----------
        self: CapturedCs loading instance
        
        str_inst: Stream instance
            Gas FP from melt-refining separations to be captured onto solid
        
        time (optional): str or float
            Time to grab Stream instance data to transfer to waste form
            (must correspond to dataframe column or keyword)
        
        **kwargs (optional):
        - verbose: if True, print details about waste loading
        - loading_fraction: float
            User-specified waste loading fraction
        """
        
        loaddata = self._get_loading_data(str_inst.id, str_inst.form)
        temp = deepcopy(str_inst)
        series = stream.get_srs(temp, time)
        to_drop = []
        for attr, srs in series:
            for nuc in srs.index:
                el = nuclidedata.determine_nuclide_info(nuc)[0]
                if el not in ('cs', 'rb'):
                    to_drop.append(nuc)
            srs = srs.drop(to_drop)
            setattr(temp, attr, srs.to_frame())
        temp.mass = sum(temp.comp[stream._get_time(time, temp.comp)]) #g cs,rb
        if(kwargs.get('loading_fraction')):
            n_can = self.easy_hlw_loading(temp, kwargs['loading_fraction'])
        else:
            n_can = temp.mass/loaddata.alkali_per_canister # expect << 1
        return(self._make_wasteform(temp, n_can, time, **kwargs))


################################################################################
# WASTE CRUCIBLE SKULL/SLAG FROM MELT REFINING
################################################################################


class Skull(AqGlass):
    
    """
    In melt refining, some of the elements in the melt form a skull on the
    crucible. After the melt is poured off to be recast into new fuel, the 
    skull is oxidized (burned?), at which point it can be removed from the
    crucible and made into a wasteform.

    If melt refining is to be used in a continuous recycle process, the skull
    would need to be processed to further recover U and Pu. However, because 
    melt refining is only used in the limited-recycle options, and its use is 
    limited to 3 recycle steps, it is assumed that no skull processing is 
    required. Therefore, I assume the skull is removed from the crucible via 
    oxidation with no postprocessing. The oxidized skull is then made into a 
    glass waste form similar to UREX glass for disposal
    
    """
    
    def __init__(self, **kwargs):
        super(Skull, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)


################################################################################


