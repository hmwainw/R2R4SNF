###############################################################################
"""
    Last edited on February 17, 2019

    @author: matz

    comments: Separation class and methods
    
"""
###############################################################################
import os
import imp
import numpy as np
import pandas as pd
from copy import deepcopy
from collections import Counter
import itertools
from nwpy import stream
from nwpy import nuclidedata
###############################################################################
# SOLID FUEL SEPARATION
###############################################################################


class Separation(object):
    
    
    """
    In most cases, the fuel cycle data reflects solid fuels, for which
    separations are calculated using the following class. The calculation
    is carried out by the method "reprocess", which takes as input the
    elements to recover as products, the recovery fractions (i.e. what 
    percent of that element in the feed is recovered), and the name/type
    of process with which the separation is performed.
    
    The calculation is split into two parts which are both contained in 
    the main method "reprocess". The first submethod, "recover", uses 
    the recovered elements and fractions as input to remove them from 
    the feed stream as products.
    
    Then, the remaining material is processed by the second method,
    "partition", which takes as input the name/type of separation process.
    With that, a data file informs the method about the waste streams it
    produces.
    
    """
    
    
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.method = self.data['reprocessing']
        self.recovered = self.data['recovered']
        
        
    #--------------------------------------------------------------------------
    # REPROCESS
    # Recover actinides from used fuel stream and return effluent waste
    # streams for further processing and waste form loading
    #--------------------------------------------------------------------------
    def reprocess(self, str_inst, rfrac=0.99, rdict={}, time='last', **kwargs):
        """Recover actinides from used fuel stream and return effluent waste
        streams for further processing and waste form loading
        
        Parameters
        ----------
        str_inst: Stream instance
        
        rfrac (optional): float
            Fractional recovery of elements indicated in stage data
        
        rdict (optional): dict
            Dictionary containing fractional recovery for specific elements;
            all others to be treated as default
            
        time (optional): str or float
            Indicates time at which to impose separation on str_inst
            
        Returns
        -------
        Waste stream instance
        
        """
        
        recovery = self._build_recovery_dict(rfrac, rdict, **kwargs)
        product_strm, waste_strm = self.recover(str_inst, recovery, time)
        try:
            wastes = self.partition(waste_strm, time)
        except: # expected msr behavior
            wastes = waste_strm
        return(wastes)
    
    
    def _build_recovery_dict(self, frac, user_sp_dict={}, **kwargs):
        """Make the dictionary that defines recovery of product species;
        inputs are pulled from fuel cycle data (recovered elements), 
        separation process data, and user input.
        
        Parameters
        ----------
        frac: float
            Recovery fraction for recovered elements
        
        user_sp_dict (optional): dict
            Dictionary input by user indicating elements and recovery
            fractions not covered by the above fraction
        
        kwargs:
        
        Returns
        -------
        Dictionary of recovery fractions for recovered elements
        
        """
        
        d = {}
        for sp in self.recovered.split(','): # first, get from .fc
            d = self._update_recovery_dict(d, sp, frac)
        if(self.method in ['meltrefining', 'msr']): # then, get from .py
            d = self._get_special_recovery_dict(d, **kwargs)
        for sp in user_sp_dict.keys(): # then, get from user
            d = self._update_recovery_dict(d, sp, user_sp_dict[sp])
        return(d)
        
    
    @staticmethod
    def _update_recovery_dict(r_dict, species, fraction):
        """Split group specifications into elements, update dictionary"""
    
        if(species in ['act', 'actinide', 'actinides']):
            to_append = nuclidedata.actinides
        elif(species in ['tru', 'transuranic', 'transuranics']):
            to_append = nuclidedata.transuranics
        elif(species in ['ma', 'minor actinide', 'minor actinides',
                    'minoractinide', 'minoractinides']):
            to_append = nuclidedata.minor_actinides
        elif(species in ['u3']):
            to_append = nuclidedata.u3
        else: # append the species str directly
            to_append = [species]
        for spi in to_append:
            r_dict[spi] = fraction
        return(r_dict)
    
    
    def _get_special_recovery_dict(self, r_dict, **kwargs):
        """Once product recovery dict has been made, account for speciation 
        of non-actinides among products and wastes in melt-refining or msr
        separations processes
        
        Parameters
        ----------
        r_dict: dict
            Dictionary indicating recovery fractions (1-x_eff) for user-
            requested elements
        
        kwargs:
            - 're': if ='lo', rare-earth element separation efficiencies are
            10x lower than the default
        
        Returns
        -------
        Updated recovery dictionary containing separation-specific recoveries
        for non-user-requested species (such as FP that are recovered by means
        of the given separations process
        
        """
        
        # overwrites previous recovery dict
        # get other dict from datapath
        if(kwargs.get('re')=='lo'):
            file = self.method+'_re-lo'+'.py'
        else:
            file = self.method+'.py'
        sep = imp.load_source(file, os.path.join(self.datapath, 'sep', file))
        new = sep.recover
        for sp in new.keys(): # overwrite r_dict with fractions from .py
            r_dict[sp] = new[sp]
        return(r_dict)
        
        
#        for sp in new.keys(): # do not overwrite requested elements
#            if(sp in self.recovered.split(',')):
#                new.pop(sp, None)
#        for species in r_dict:
#            if(species in new):
#                continue
#            else:
#                new[species] = r_dict[species]
#        return(new)


    # This method is used in both the recover and
    # partition methods so it is defined up here
    @staticmethod
    def split_series(srs, sep):
        """Split a Pandas Series object into one or more new
        Series according to specified fractions
        
        Parameters
        ----------
        srs: Pandas Series object
            The Series to be split
    
        sep: NumPy Array
            The array that controls the splitting of the Series
            The rows represent elements; the cols represent product streams
            
        Returns
        -------
        List of new Pandas Series
        
        """
        
        arr = np.zeros((len(srs), sep.shape[1]))
        for i in range(0, len(srs)):
            nuc = srs.index[i]
            z = nuclidedata.Z[nuclidedata.determine_nuclide_info(nuc)[0]]-1
            arr[i] = srs[nuc]*sep[z]
        new_srs = [pd.Series(arr[:,j], index=srs.index, name=srs.name)
                   for j in range(0, arr.shape[1])]
        return(new_srs)


    #--------------------------------------------------------------------------
    # RECOVER
    # In this method, the recovered species are removed from the feed stream.
    # The method iterates through the feed dataframe, partitioning isotopes
    # into product and wastes. These are kept as dictionaries and lists. Then
    # new dataframes are made and new Streams produced
    #--------------------------------------------------------------------------
    def recover(self, str_inst, r_dict, time):
        """Recover species from feed stream to make product stream
        
        Parameters
        ----------
        self: Separation instance
        
        str_inst: Stream instance
            The feed to be reprocessed
        
        r_dict: dict
            Keys correspond to the species (str)
            Values to the recovery fraction (float)
            
        time:
        
        Returns
        -------
        Separated product and waste Stream objects
        
        """
        
        # preserve other attributes such as eg, stgn, etc in outstreams
        # in outstreams list, item 1 is product stream; item 2 is waste stream
        outstreams = [deepcopy(str_inst), deepcopy(str_inst)]
        series = stream.get_srs(str_inst, time)
        sep_eff = self._make_separation_array(r_dict)
        for attr, old_srs in series:
            new_srs_list = self.split_series(old_srs, sep_eff)
            for i in range(0, len(new_srs_list)):
                setattr(outstreams[i], attr, new_srs_list[i].to_frame())
                if(attr=='comp'):
                    outstreams[i].mass = sum(new_srs_list[i])
        return(self._finalize_recovered_streams(str_inst, outstreams))



    @staticmethod
    def _finalize_recovered_streams(instream, outstrm):
        """Modify some of the attributes of the product and waste
        streams produced after recovering actinides"""
        
        prod, waste = outstrm
        del prod.id
        prod.form = 'products'
        waste.id = 'hlw'
        waste.form = 'wastes'
        if(any([x for x in ['driver', 'blanket'] if x==instream.form])):
            waste._dfbb = instream.form # hidden attribute id
        return([prod, waste])
    
    
    @staticmethod
    def _make_separation_array(rd):
        """With the input dictionary of elements and recovery fractions, 
        make an array to apply to the separation of streams"""
        
        # rows: elements, cols: (products, wastes)
        # elements first; element specification > group specification
        sep = np.concatenate((np.zeros((99,1)), np.ones((99,1))), axis=1)
        hist = []
        for i in range(0, len(sep)):
            el = list(nuclidedata.Z.keys())[list(nuclidedata.Z.values()).index(i+1)]
            matches = [nuclidedata.is_group(el, k) for k in rd.keys()]
            if(any(matches)):
                k = list(rd.keys())[np.where(matches)[0][0]]
                sep[i] = np.array((rd[k], 1-rd[k]))
        return(sep)


    #--------------------------------------------------------------------------
    # PARTITION
    # In this method, the waste stream from "recover" is partitioned into
    # multiple waste streams reflecting expected effluents from the specified
    # process. The characteristics (mass, composition, decay heat, etc) are
    # determined for each effluent stream.
    #--------------------------------------------------------------------------
    def partition(self, str_inst, time):
        """Get data file with separation efficiencies and apply to
        waste stream, returning new waste streams and compositions
        
        Parameters
        ----------
        str_inst: Stream instance
            The wastes from reprocessing the feed stream
            
        Returns
        -------
        Waste Stream instances resulting from reprocessing
        
        """
        
        sep, strm_names = self._read_sep_data()
        outstreams = [stream.empty() for i in range(0, sep.shape[1])]
        series = stream.get_srs(str_inst, time)
        for name, old_srs in series:
            new_srs_list = self.split_series(old_srs, sep)
            for i in range(0, len(new_srs_list)):
                setattr(outstreams[i], name, new_srs_list[i].to_frame())
        out = self._finalize_waste_streams(str_inst, outstreams, strm_names)
        return(out)


    @staticmethod
    def _finalize_waste_streams(instream, outstreams, strm_names):
        """Modify some of the attributes of the waste
        streams produced after partitioning"""
        
        for j in range(0, len(outstreams)):
            outstreams[j].form = strm_names[j].lower()
            if(hasattr(instream, '_dfbb')):
                # append "driver" or "blanket" to clarify waste forms
                outstreams[j].form += '_'+instream._dfbb
            outstreams[j].id = 'hlw'
            t = outstreams[j].comp.columns[0] # only one col
            outstreams[j].mass = sum(outstreams[j].comp[t])
        # misc additional attributes
        for a in ['batch', 'evaluationgroup', 'stagenumber']:
            if(hasattr(instream, a)):
                for j in range(0, len(outstreams)):
                    setattr(outstreams[j], a, getattr(instream, a))
        return(outstreams)
        
        
    def _read_sep_data(self):
        """Import separation efficiency data from the sep directory"""
            
        sep_data, path = self._import_separation_data()
        outputs, hlines = self._get_output_streams(sep_data)
        cols = tuple(np.arange(2, len(outputs)+2))
        sep_data = np.loadtxt(path, skiprows=hlines, usecols=cols)
        names = [outputs[i-2] for i in cols]
        return(sep_data, names)


    def _import_separation_data(self):
        """Open data file containing evaluation group separation data"""
            
        file = self.method+'.sep'
        dp = os.path.join(self.datapath, 'sep', file)
        sep_data = open(dp).read().splitlines()
        sep_data = filter(None, sep_data)
        # sep_data = itertools.filterfalse(None, sep_data)
        sep_data = [x for x in sep_data]
        return(sep_data, dp)


    @staticmethod
    def _get_output_streams(sep_data):
        """Determine waste forms, count datafile header lines"""
        
        header=[]
        header_lines = 0
        for line in sep_data:
            if(int(line.split()[0])==0):
                header.append(line)
                header_lines += 1
        # stream names on last line of header
        return(header[-1].split()[2:], header_lines)


###############################################################################
# LIQUID FUEL SEPARATION
###############################################################################
#For MSR
#- discharge, cool, apply separation
#- scale to reflect mass for one year of operation
#- run origen calculation reporting values for each year of operation
#- gives results at [discharge, sep, 1, 2, 3, ..., n-2, n-1, n] years
#    These are the wastes accumulated by MSR
#        Total canisters, canisters per year
#        All loaded the same, just with different decay times

class LiquidFuelSep(object):
    
    
    """
    In some cases (namely, MSR), the data reflect the use of liquid fuel that
    flows continuously in and out of the core - in those cases, this function
    is used to apply separations to calculate the output material streams.
    
    
    
    """
    
    
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.wastespecies = self.data['wastespecies'].lower().split(',')
        self.dischargedfuel = self.data['dischargedfuel']
        self.wastemass = self.data['wastemass']
        

    def reprocess(self, streams):
        """For liquid (MSR) fuel, apply separation and account
        for unrecovered actinides
        
        
        Parameters
        ----------
        streams: list
            List containing two streams instances
            1. Salt discharged from the MSR and diverted to MSR separations
            2. Fission product stream separated from salt in MSR separations
        
        Returns
        -------
        Stream instances resulting from reprocessing
        
        """
        
        outstreams = []
        temp = deepcopy(streams)
        df_str, fp_str = temp
        assert df_str.form=='df1', 'DF Stream form given as '+str(df_str.form)
        assert fp_str.form=='fp2', 'FP Stream form given as '+str(fp_str.form)
        # fission product stream
        fp_str = self._get_fp_stream_mass(fp_str)
        fp_str = self._makeup(fp_str, df_str)
        fp_str = self._apply_tol(fp_str)
        # discharged salt stream (if applicable)
        if(self.data['dischargedfuel']!=0.0):
            df_str = self._apply_tol(df_str)
            df_str.comp = df_str.comp*self.dischargedfuel*1e6/df_str.mass
            df_str.heat = df_str.heat*self.dischargedfuel*1e6/df_str.mass
            outstreams.append(Stream(self.dischargedfuel*1e6, comp=df_str.comp,
                                     heat=df_str.heat, time=df_str.time,
                                     evaluationgroup=df_str.evaluationgroup,
                                     stagenumber=df_str.stagenumber,
                                     id='snf', form='dischargedsalt'))
        outstreams.append(fp_str)
        return(outstreams)


    def _get_fp_stream_mass(self, fp_stream):
        """Calculate the mass of the fp stream from msr separations"""
        
        # calculate FP fraction
        x_fp = fp_stream.mass_fraction('fp')
        # scale fp mass, composition by expected mass of fp
        m_fp = self.wastemass[self.wastespecies.index('fp')] # tons
        fp_stream.comp = fp_stream.comp*m_fp/x_fp
        fp_stream.heat = fp_stream.heat*m_fp/x_fp
        fp_stream.mass = fp_stream.comp[fp_stream.comp.columns[-1]].sum()
        return(fp_stream)


    def _makeup(self, w_stream, df_stream):
        """For a given species in the waste, determine if isotopic composition
        accounts for its expected mass; if not, make up using SNF stream"""
        
        for i in range(0, len(self.wastespecies)):
            species = self.wastespecies[i]
            mass_in_w = w_stream.mass*w_stream.mass_fraction(species)
            diff = self.wastemass[i]*1e6 - mass_in_w
            if(diff > 0.0):
                mass_in_df = df_stream.mass*df_stream.mass_fraction(species)
                multiplier = diff/mass_in_df
                for nuc in df_stream.comp.index:
                    if(nuclidedata.is_group(nuc, species)):
                        if(nuc in w_stream.comp.keys()):
                            w_stream.comp[nuc]+=df_stream.comp[nuc]*multiplier
                        else:
                            w_stream.comp[nuc]=df_stream.comp[nuc]*multiplier
                        if(nuc in w_stream.heat.keys()):
                            w_stream.heat[nuc]+=df_stream.heat[nuc]*multiplier
                        else:
                            w_stream.heat[nuc]=df_stream.heat[nuc]*multiplier
        w_stream.mass = sum(w_stream.comp.values())
        return(w_stream)


    def _apply_tol(self, stream):
        """Apply the mass cutoff to a stream"""
        
        to_delete = []
        for nuc in stream.comp:
            if(stream.comp[nuc] < self.tol):
                to_delete.append(nuc)
        for ignored in to_delete:
            del stream.comp[ignored]
            try:
                del stream.heat[ignored]
            except:
                pass
        stream.mass = sum(stream.comp.values())
        return(stream)


###############################################################################

