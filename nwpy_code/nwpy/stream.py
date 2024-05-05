###############################################################################
"""
    Last edited on May 6, 2019

    @author: matz

    comments: Stream and WasteForm instances for fuel cycle calculations
    
"""
###############################################################################
import os
import pandas as pd
from nwpy import nuclidedata
#from plot import TimePlot as tp
#from plot import PiePlot as pp
###############################################################################
# STREAM CLASS
###############################################################################


class Stream(object):
    
    
    """
    Instances of the Stream class represent the flows of material between 
    fuel cycle stage operations such as irradiation systems, cooling, and
    reprocessing.
    
    Fundamentally, the Stream is only a container for objects. The objects
    must be named a certain way for other methods in fuelcycle to find them.
    On instantiation, the only requirements for the Stream are that it have 
    a non-negative mass.
    
    """
    
    
    # specify the units of the properties utilized in ORIGEN calculations
    units = {}
    units['grams'] = 'comp'     # composition
    units['watts'] = 'heat'     # decay heat
    units['Ci'] = 'act'         # radioactivity
    units['m3 water'] = 'tox'   # radiotoxicity
    
    
    def __init__(self, mass, **kwargs):
        """Initialize the Stream instance, test the input values"""
        
        assert isinstance(mass, (int, float)), 'Mass must be a number'
        assert mass >= 0.0, 'Stream mass must be non-negative.'
        self.mass = mass
        self.form = None # will be overwritten if in kwargs
        for key, value in kwargs.items():
            # expected keys:
            # comp, heat, act, tox (pd.DataFrame)
            # time (float), form (str), batch (int), id (int)
            setattr(self, key, value)



    def __repr__(self):
        return('Stream instance: '+str(self.form))
        
        
    def __str__(self):
        p = str(self.mass)+' grams '+str(self.form)
        return(p)
    
    
    #--------------------------------------------------------------------------
    # ADD AND SUBTRACT STREAMS
    # With 'add', two streams are combined, summing their masses and the
    # values for all associated data (such as composition, heat, etc.) With
    # 'subtract', one stream is deducted from another.
    #
    # Note: The 'subtract' method may be unstable and should be used
    # carefully. Although one stream may have greater mass than the other,
    # negative values can result in the data attributes (for composition,
    # heat, etc). One way to mitigate this risk is to only subtract streams
    # you know are related.
    #--------------------------------------------------------------------------
    def add(self, other, time='last', id=None, form=None):
        """ Combine two streams using the composition at a specific time 
        
        Parameters
        ----------
        self: Stream object
        
        other: Stream object
            The Stream to add
        
        time: float or str
            Time (float) to pull and sum composition data from both streams
            If 'last' (str), use the data from the last column entry.
            If 'first' (str), use data from the first column entry.
        Returns
        -------
        Stream object
        
        """
        
        if(any([x.form == 'empty' for x in [self, other]])): # return nonempty
            return([self, other][1-[self.form, other.form].index('empty')])
        if(isinstance(time, str)):
            time = self._parse_time(time)
        self._check_time(other, time)
        if(id == None and hasattr(self, 'id')): # maintain id from self
            id = self.id
        if(form == None and hasattr(self, 'form')): # maintain form from self
            form = self.form
        new = Stream(self.mass+other.mass, id=id, form=form)
        # combine dataframes
        for x in dir(self):
            if(x in dir(other) and all([isinstance(getattr(j, x), pd.DataFrame)
                                        for j in [self, other]])):
                df = [getattr(self, x), getattr(other, x)]
                #if(any([i.empty for i in df])):
                #    idx=[i for i,j in enumerate(df) if j.empty][0]
                #    df[idx] = pd.DataFrame({time: []})
                try:
                    new_df = df[0][time].add(df[1][time], fill_value=0.0)
                except:
                    print(x)
                    raise
                new_df = new_df.to_frame(name=time)
                setattr(new, x, new_df)
        # evaluation group
        if(self._compare_eg(self, other)):
            new.evaluationgroup = self.evaluationgroup
        # stage number
        if(self._compare_stgn(self, other)):
            new.stagenumber = self.stagenumber
        # return the new stream
        return(new)


    def subtract(self, other, time='last', id=None, form=None):
        """Subtract one Stream object from another at a specified time.
        
        Parameters
        ----------
        
        self: Stream object
            
        other: Stream object
            This is the Stream that is subtracted from self
            
        time: float
            This is the time at which the composition data is returned; must
            reflect a value in the columns of the Stream dataframes.
        
        Returns
        -------
        Stream object
        
        """
        
        if(any([x.form == 'empty' for x in [self, other]])): # return nonempty
            return([self, other][1-[self.form, other.form].index('empty')])
        assert self.mass > other.mass, 'Stream1 - Stream2 has negative mass'
        if(isinstance(time, str)):
            time = self._parse_time(time)
        self._check_time(other, time)
        if(id==None): # maintain id from self
            id = self.id
        if(form==None): # maintain form from self
            form = self.form
        new = Stream(self.mass-other.mass, id=id, form=form)
        # combine dataframes
        for x in dir(self):
            if(x in dir(other) and all([isinstance(getattr(j, x),pd.DataFrame)
                                        for j in [self, other]])):
                df = [getattr(self, x), getattr(other, x)]
                #if(any([i.empty for i in df])): # assume only one is empty
                #    idx=[i for i,j in enumerate(df) if j.empty][0]
                #    df[idx] = pd.DataFrame({time: []})
                new_df = df[0][time].subtract(df[1][time], fill_value=0.0)
                new_df = self._check_negatives(new_df, x)
                new_df = new_df.to_frame(name=time)
                setattr(new, x, new_df)
        # evaluation group
        if(self._compare_eg(self, other)):
            new.evaluationgroup = self.evaluationgroup
        # stage number
        if(self._compare_stgn(self, other)):
            new.stagenumber = self.stagenumber
        # return the new stream
        return(new)


    def _parse_time(self, t_str):
        """If time is a string, figure out how to process it"""
        
        if(t_str == 'last'):
            t_str = self.comp.columns[-1]
        elif(t_str == 'first'):
            t_str = self.comp.columns[0]
        else:
            try:
                t_str = float(t_str)
            except:
                raise
        return(t_str)


    def _check_time(self, other, t):
        """Ensure that specified time appears in data for both streams"""
    
        for x in dir(self):
            if(x in dir(other) and all([isinstance(getattr(j, x),pd.DataFrame)
                                        for j in [self, other]])):
                temp1 = getattr(self, x)
                temp2 = getattr(other, x)
                if(any(df.empty for df in [temp1, temp2])):
                    return
                assert all(t in c for c in [temp1.columns,temp2.columns]),\
                    'Specified time not in attribute '+x


    @staticmethod
    def _check_negatives(srs, var='mass'):
        """Check for negative values in Pandas Series; warn user, set to 0
        
        Parameters
        ----------
        srs: Pandas Series
            Contains composition, heat, etc. data
            
        var (optional): str
            Name of variable, e.g. 'heat'
        
        Returns
        -------
        Non-negative Pandas Series
        
        """

        for j in srs.index:
            if(srs[j] < 0.0):
                print('Negative '+var+' for nuclide '+j+'; setting to 0.0')
                srs[j] = 0.0
        return(srs)


    @staticmethod
    def _compare_eg(s1, s2):
        """Compare the evaluationgroup attribute for two Streams"""
        
        if(all([hasattr(x, 'evaluationgroup') for x in [s1, s2]])):
            if(s1.evaluationgroup == s2.evaluationgroup):
                return(True)
        return(False)


    @staticmethod
    def _compare_stgn(s1, s2):
        """Compare the stagenumber attribute for two Streams"""
        
        if(all([hasattr(x, 'stagenumber') for x in [s1, s2]])):
            if(s1.stagenumber == s2.stagenumber):
                return True
        return False
    
    
    #--------------------------------------------------------------------------
    # UPDATE
    # Append new columns to one of the Pandas DataFrame attributes and sort by
    # columns (times)
    #--------------------------------------------------------------------------
    def update(self, attr, df, time_to_add='last', overwrite=False):
        """Add or update Stream DataFrame attributes and sort by columns
        
        Parameters
        ----------
        self: Stream object
        
        attr: str 
            Attribute of Stream object to be created or updated
        
        df: DataFrame
            Pandas DataFrame to add as attribute or append to existing 
            attribute
        
        time_to_add (optional): float
            ORIGEN-S returns values at times starting from 0.0; if the
            input was made at a later time, the times must be shifted to
            reflect that. The default value is a string, 'last', which 
            indicates the last time in the DF for which there is data. If
            the Stream does not have the requested attribute, the value is 0.0
            
        Returns
        -------
        None
        
        """
        
        if(not hasattr(self, attr) or overwrite==True):
            setattr(self, attr, df)
        else:
            dat = getattr(self, attr)
            if(time_to_add=='last'): # shift new df columns
                time_to_add = dat.columns[-1]
            df.columns = df.columns+time_to_add
            try:
                dat = dat.join(df, how='outer')
            except ValueError: # duplicate columns
                setattr(self, attr, df)
            else:
                dat = dat.fillna(0.0)
                dat = dat.sort_index(axis=1)
                setattr(self, attr, dat)


    #--------------------------------------------------------------------------
    # MASS FRACTION CALCULATOR
    # Calculate the mass fraction of isotopes, elements, or groups of elements
    # (i.e. actinides) with this method.
    #--------------------------------------------------------------------------
    def mass_fraction(self, species, time='last'):
        """Return the mass fraction of species in a stream
            
        Parameters
        ----------
        self: Stream instance
        
        species: str
            Isotope or element or species ('FP', 'actinide', 'TRU', etc)
            
        Returns
        -------
        float: Mass fraction of requested species
        
        """
        
        if(time=='last'):
            time = self.comp.columns[-1]
        species = species.lower()
        mf = 0.0
        try: # try species as isotope
            mf += self.comp[time][species]
        except KeyError:
            # can only do it this way because these methods have no overlap
            for nuclide in self.comp.index: # iterate over pandas df index
                if(nuclidedata.is_group(nuclide, species)):
                    mf += self.comp[time][nuclide] # treat species as group
                elif(nuclidedata.determine_nuclide_info(nuclide)[0]==species):
                    mf += self.comp[time][nuclide] # treat species as element
                else:
                    continue
        return(mf/self.mass)


    #--------------------------------------------------------------------------
    # WRITING RESULTS TO FILE
    # Make the results available to other codes, especially those running
    # remotely, by writing them to a file that can be transplanted
    #--------------------------------------------------------------------------
    def write(self, path, property='comp', total=True, **kwargs):
        """Write the stream data and its properties to a text file that
        can be accessed by other codes. Currently only writes the stream
        total for time-dependent data; in the future, will write composition-
        dependent data.
        
        Parameters
        ----------
        self: Stream instance
            The object containing the waste to decay
    
        path: str
            Indicate a directory path to write the file
        
        property (optional): str
            Stream attribute to write to file
            - comp (composition)
            - heat
            - act (radioactivity)
            - tox (radiotoxicity)
        
        total (optional): Boolean
            If True, the data to be written is the stream total 
            (summing the data for all isotopes for each time point)
            
        kwargs
        ------
        name: str
            File name
        evaluationgroup: str
            Fuel cycle evaluation group (to be included in data file)
        stagenumber: str
            Fuel cycle stage number (to be included in data file)
        
        Returns
        -------
        Path to the file that was written
        
        """
        
        filename = self._make_filename(**kwargs)
        if(total is not True):
            raise NotImplementedError('tbd')
        dat = ''
        if(kwargs.get('name')):
            dat += 'name,'+kwargs['name']+'\n'
        if(kwargs.get('evaluationgroup')):
            dat += 'evaluationgroup,'+kwargs['evaluationgroup']+'\n'
        if(kwargs.get('stagenumber')):
            dat += 'stagenumber,'+kwargs['stagenumber']+'\n'
        dat += 'id,'+self.id+'\n'
        dat += 'form,'+self.form + '\n'
        if(hasattr(self, 'canister')):
            dat += 'pkg_diameter,'+str(self.canister['Diameter'])+'\n'
            dat += 'pkg_length,'+str(self.canister['Length'])+'\n'
        if(hasattr(self, 'number')):
            dat += 'n_wf,'+str(self.number)+ '\n'
        # Data: sum temp for all times
        temp = getattr(self, property)
        temp = temp.sum(axis=0)
        for t in temp.index:
            dat += str(t)+','+str(temp[t])+'\n'
        file = open(os.path.join(path, filename), 'w')
        file.write(dat)
        file.close()
        # report the total mass? total canisters?
        # wasteforms per canister?


    def _make_filename(self, **kwargs):
        """Make the data file name"""
        
        filename = ''
        if(kwargs.get('name')):
            return(kwargs['name'])
        if(kwargs.get('evaluationgroup')):
            filename += kwargs['evaluationgroup']
        if(kwargs.get('stagenumber')):
            filename += kwargs['stagenumber']
        filename += self.form
        return(filename+'.csv')


###############################################################################
# MODULE METHODS
###############################################################################


#------------------------------------------------------------------------------
# EMPTY STREAM
# Make an empty Stream object; one example of potential use is as a
# container to which other Streams can be added.
#------------------------------------------------------------------------------
def empty():
    """Make an empty Stream container"""
    
    return(Stream(0.0, form='empty'))


#------------------------------------------------------------------------------
# GET DATAFRAMES
# Get a list of all dataframes and the attribute names contained in a Stream
# object; returned as tuples with name, df pairs
#------------------------------------------------------------------------------
def get_df(strm):
    """Get the attribute name and value for all DataFrame objects in a Stream"""
    
    return([(i, j) for i,j in strm.__dict__.items()
            if isinstance(j, pd.DataFrame)])


#------------------------------------------------------------------------------
# GET SERIES
# Get a list of all series and the attribute names contained in a Stream
# object; returned as tuples with name, series pairs; the Series is the subset
# of the DataFrame by the column header
#------------------------------------------------------------------------------
def get_srs(strm, time): # get srs nerd
    """Return name and Series subset by column header (time) for all 
    DataFrame objects in a Stream
    
    Parameters
    ----------
    strm: Stream instance
    
    time: kw
        Time argument to subset DataFrames. Options include
        - any float or int
        - 'last' (str): return the last column
        - 'first' (str): return the first column
    
    Returns
    -------
    List of tuples containing (attribute name, Pandas Series object)
    
    """

    df_list = get_df(strm)
    srs_list = []
    for name, df in df_list:
        t = _get_time(time, df)
        srs_list.append((name, df[t]))
    return(srs_list)


def _get_time(t, df):
    """Take different keyword inputs for time and return the appropriate 
    column name from the Pandas DF given."""
    
    if(t not in df.columns):
        if(t == 'last'):
            t = df.columns[-1]
        elif(t == 'first'):
            t = df.columns[0]
        else:
            try:
                t = float(t)
            except:
                raise
    assert(t in df.columns)
    return(t)


###############################################################################
# WASTE CLASS
###############################################################################


class WasteForm(Stream):
    
    
    """
    The WasteForm class is a subclass of the Stream class with some extra
    attributes. WasteForm instances are produced within the Stage class in 
    the load_waste method.
    
    """
    
    def __init__(self, mass, number, canister, **kwargs):
        """Initialize the WasteForm instance; the WasteForm requires the
        same input as the Stream (only the mass) and also information about
        the number and dimensions of the canisters."""
        
        super(WasteForm, self).__init__(mass, **kwargs)
        self.number = number    # number of canisters
        self.canister = canister # dict; dimensions

    
    def __repr__(self):
        return('Wasteform instance: '+str(self.form))
    
    
    def __str__(self):
        p = str(self.form)+' waste form canister'+'\n'
        p += str(round(self.mass/1e3, 2))+' kg waste per canister'+'\n'
        #p += str(round(sum(self.heat.values()), 2))+ ' W per canister'+'\n'
        try: # snf canisters don't have this
            mass_limit = round(self.canister['Mass limit']/1e3, 2)
        except:
            pass
        else:
            p += ('Canister mass limit (kg; waste + matrix): '+
                  str(mass_limit)+'\n')
        p += str(self.number)+' canisters'
        return(p)


###############################################################################