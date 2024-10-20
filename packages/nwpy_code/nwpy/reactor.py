###############################################################################
"""
    Last edited on May 13, 2019

    @author: matz

    comments: Contains the Reactor object as well as special case classes
    
"""
###############################################################################
import os
import pandas as pd
from nwpy.stream import Stream
from nwpy import nuclidedata
###############################################################################


class Reactor(object):
    
    
    """
    The Reactor object returns to the Stage instance the isotopic composition
    of the stream discharged from the irradiation system in that stage. The
    format of the isotopic csv files is based on the type of irradiation
    system. The base class for Reactor objects; works for all general cases, 
    including PWR, SFR (breed or burn), HTGR, HWR, or EDS.
    
    The methods in the Reactor class require many instantiated attributes to
    be passed in from the Stage instance. In particular, the stage data and
    the datapath are required in order to properly find and handle the
    composition data.
    
    """
    
    
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


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
        
        c = self.get_composition_data()
        outstreams = [] # list to hold output streams
        for strm in c: # mass in grams; here, strm is stream name
            m = sum(c[strm])
            df = c[strm].to_frame(name = self.get_time())
            temp_strm = Stream(m, comp=df, form=strm, id='snf',
                               evaluationgroup=self.evaluationgroup,
                               stagenumber=self.number)
            if(not any([True for x in ['snf', 'batch', 'df',
                                       'driver','blanket'] if x in strm])):
                temp_strm.id = 'hlw'
            if('batch' in strm):
                try:
                    temp_strm.batch = int(strm[5:])
                except ValueError:
                    temp_strm.batch = int(strm[5:-3]) # remove "_df"
            outstreams.append(temp_strm)
        return(self._make_outstreams(outstreams))
    
    
    def get_composition_data(self):
        """Depending on the reactor type, get the isotopic data
    
        Parameters
        ----------
        self: Reactor instance
    
        Results
        -------
        Dictionary of Pandas Series, which contain the composition of 
            each stream produced in the Stage
    
        """
        
        comp = {}
        df = self._import_isotopic_csv()
        cols, masses, labels = self._streams(df=df) # mass in t
        df = df[cols]
        for i in range(0, len(df.columns)):
            column = df.columns[i]
            renorm = 1e6/sum(df[column])
            label = labels[i]
            if(masses[i]<0.0):
                label = label+'_df' # for breed and burn sfr batches
            comp[label] = self._rescale_df(df[column], abs(masses[i]), renorm)
            comp[label] = nuclidedata.group_daughters(comp[label])
        return(comp)


    def get_time(self):
        """Specify the time at which the stream was discharged"""
    
        if(self.data['coolingtime']!=0.0):
            t = 0.0
        else:
            t = 5.0
        return(t)
    
    
    @staticmethod
    def _make_outstreams(outstream_list):
        """Return a list of Stream instances or a single Stream instance"""
        
        if(len(outstream_list)==1):
            return(outstream_list[0])
        else:
            return(outstream_list)


    def _streams(self, **kwargs):
        """Info to read in isotopic data for general evaluation groups"""
        
        df = kwargs['df']
        stream_masses = [self.data['mass']]
        labels = ['snf']
        col_names = list(df.columns[df.columns.str.contains('discharge')])
        col_names = [c for c in col_names if str('ST'+str(self.number)) in c]
        return(col_names, stream_masses, labels)
    
    
    @staticmethod
    def _rescale_df(srs, mass, x=1.0):
        """Rescale the values in the Pandas Series to sum to the stream mass"""
        
        # reindex the dataframe so the strs are lowercase
        srs.index = srs.index.str.lower()
        # rescale the values to sum to stream mass
        srs.loc[:] *= mass*x
        return(srs)
        
    
#    @staticmethod
#    def _group_daughters(srs):
#        """Group nuclides unsupported by ORIGEN with their decay daughters"""
#        
#        to_drop = []
#        for j in range(0, len(srs)):
#            nuc = srs.index[j].lower()
#            mass = srs[srs.index[j]]
#            if(nuc in nuclidedata.group_nuclides.keys()):
#                to_drop.append(nuc)
#                for daughter in nuclidedata.group_nuclides[nuc].keys():
#                    branchfrac = nuclidedata.group_nuclides[nuc][daughter]
#                    try: # assume daughter is already in composition df
#                        srs[daughter] += mass*branchfrac
#                    except KeyError: # if daughter not in df, make new entry
#                        srs[daughter] = mass*branchfrac
#        srs = srs.drop(to_drop)
#        return(srs)


    def _import_isotopic_csv(self):
        """Open csv file for the isotopic data for the evaluation group"""
        
        file = os.path.join(self.datapath, 'iso', self.evaluationgroup+'.csv')
        return(pd.read_csv(file, index_col='isotope'))


###############################################################################
# FUSION FISSION HYBRID
###############################################################################


class DFBB(Reactor):


    """
    Special case Reactor object for when masses and compositions of driver
    fuel and breeding blanket are specified separately.

    """

    def __init__(self, **kwargs):
        super(DFBB, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)


    def _streams(self, **kwargs):
        """Info to read in isotopic data for df/bb evaluation groups"""
            
        df = kwargs['df']
        stream_masses = []
        labels = []
        col_names = list(df.columns[df.columns.str.contains('discharge')])
        col_names = [c for c in col_names if str('ST'+str(self.number)) in c]
        for column in col_names:
            if('-BB-' in column):
                stream_masses.append(self.data['blanket'])
                labels.append('blanket')
            elif('-DF-' in column):
                stream_masses.append(self.data['driver'])
                labels.append('driver')
            else:
                continue
        return(col_names, stream_masses, labels)


###############################################################################
# FUSION FISSION HYBRID
###############################################################################


class FFH(Reactor):


    """
    Special case Reactor object for the molten-salt Fusion-Fission Hybrid
    reactor.
    
    """
    
    
    def __init__(self, **kwargs):
        super(FFH, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)


    def _streams(self, **kwargs):
        """Information to read isotopic data for FFH evaluation groups"""

        col_names = ['SNF (g)', 'FP1 (g)']
        labels = ['snf', 'fp1']
        stream_masses = []
        stream_masses.append(self.data['dischargedfuel'])
        stream_masses.append(self.data['salttreatment'])
        return(col_names, stream_masses, labels)


###############################################################################
# MOLTEN SALT REACTOR (MSR)
###############################################################################


class MSR(Reactor):
    
    
    """
    Special case Reactor object for liquid-fuel molten salt reactor
    
    """
    
    
    def __init__(self, **kwargs):
        super(MSR, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)


    def _streams(self, **kwargs):
        """Information to read isotopic data for MSR evaluation groups"""
        
        col_names = ['DF1 (g)', 'FP1 (g)']
        m_df1 = (self.data['mass']-self.data['salttreatment']-
                 self.data['dischargedfuel'])
        masses = [m_df1, self.data['salttreatment']]
        labels = ['fuelsalt', 'fp1']
        if(self.data['dischargedfuel'] != 0.0):
            col_names.append('SNF (g)')
            masses.append(self.data['dischargedfuel'])
            labels.append('df1') # spent fuel
        return(col_names, masses, labels)


###############################################################################
# BREED AND BURN SFR (SUSTAINABLE-SFR)
###############################################################################


class BnBSFR(Reactor):
    
    
    """
    Special case Reactor object for the breed and burn sodium fast reactor
    (otherwise known as the sustainable sodium fast reactor, SSFR)
    
    """
    
    
    def __init__(self, **kwargs):
        super(BnBSFR, self).__init__(**kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)
    

    def _streams(self, **kwargs):
        """Info to read isotopic data for B&B SFR evaluation groups"""
        
        df = kwargs['df']
        stream_masses = []
        labels = []
        col_names = list(df.columns[df.columns.str.contains('discharge')])
        for i in range(0, len(col_names)):
            labels.append('batch'+col_names[i][5])
            stream_masses.append(self.data['masses'][i])
        #    try:
        #        stream_masses.append(self.data['recycled'][i])
        #    except IndexError: # last batch would not be recycled
        #        stream_masses.append(self.data['discharged'])
        return(col_names, stream_masses, labels)


    def _make_outstreams(self, outstream_list):
        """Return a special list of Stream instances: group batches that 
        undergo separation separately from batch that is discharged"""
            
        recycled_sublist = []
        discharged_sublist = []
        for stream in outstream_list:
            stream.batch = int(stream.form[5]) # batchX
            if('_df' in stream.form): # discharged as snf
                discharged_sublist.append(stream)
                #idx = outstream_list.index(stream)
            else:
                recycled_sublist.append(stream)
        return([recycled_sublist, discharged_sublist])


###############################################################################