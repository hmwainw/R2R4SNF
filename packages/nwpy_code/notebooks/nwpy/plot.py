###############################################################################
"""
    Last edited on February 5, 2019
    
    @author: matz

    comments: Stream plotting methods
    
"""
###############################################################################
import os
import datetime
import numpy as np
import matplotlib
from matplotlib.cbook import MatplotlibDeprecationWarning
import warnings
warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
from math import atan2,degrees
import pandas as pd
import matplotlib.pyplot as plt
from decimal import Decimal
from nwpy import nuclidedata
###############################################################################
# shorthand notation for stream properties
shortcuts = {}
shortcuts['m'] = 'comp'
shortcuts['mass'] = 'comp'
shortcuts['composition'] = 'comp'
shortcuts['h'] = 'heat'
shortcuts['power'] = 'heat'
shortcuts['activity'] = 'act'
shortcuts['radioactivity'] = 'act'
shortcuts['radiotoxicity'] = 'tox'
shortcuts['toxicity'] = 'tox'
###############################################################################
# default y-axis labels for stream attributes
ylabel = {}
ylabel['comp'] = 'Mass'
ylabel['heat'] = 'Heat'
ylabel['act'] = 'Radioactivity'
ylabel['tox'] = 'Ingestion Radiotoxicity'
###############################################################################
# default y-axis units for stream attributes
unit = {}
unit['comp'] = 'g'
unit['heat'] = 'W'
unit['act'] = 'Ci'
unit['tox'] = r'm$^3$ water'
###############################################################################


class TimePlot(object):


    """
    This class builds upon the matplotlib.pyplot functionality for fuel
    cycle plotting. In particular, it enables the plotting of stream data 
    over time, which can be done after decay calculations have been made
    using the stage.decay() method. The decay calculations load Pandas 
    DataFrames into the stream object(s). Depending on what is requested 
    by the user, the code can obtain data to plot in different ways. The 
    user can use the methods of this class to compare streams and fuel
    cycles.
    
    Because the methods in this class alter the matplotlib.pylot internal
    state, it can be used in conjunction with matplotlib.pylot. If additional
    methods beyond those provided in the TimePlot class are required, a user 
    can also import matplotlib.pyplot (referred to as 'plt' from here on)
    and, after initiating a plot using TimePlot.plot, make calls to the 
    methods in plt (i.e. plt.xlim, plt.show(), plt.save(), etc).
    
    Some notes on DataFrames (fundies, if you will):
    access rows: df.loc['row name']
    access cells: df.loc['row name']['col name']
    access cols: df.loc[:]['col name']

    """
    
    streams = []
    labels = []
    property = None


    #--------------------------------------------------------------------------
    # MAIN PLOT FUNCTION
    # This function initializes the axes upon which all future modifications
    # and adjustments are performed. This function takes a list of streams and
    # plots the desired property over time for each of them.
    #--------------------------------------------------------------------------
    @classmethod
    def plot(cls, strm, property, **kwargs):
        """Plot property data against time for a Stream instance
            
        Parameters
        ----------
        cls: TimePlot class
            Class containing information about streams, plot internals, etc
        
        strm: Stream instance
            Stream or WasteForm instance containing time-dependent data for
            the property to be plotted
        
        property: str
            String indicating stream instance Pandas DF attribute containing 
            time-dependent data
            - 'comp'
            - 'heat'
            - 'radioactivity'
            - 'radiotoxicity'
        
        kwargs: Keyword formatting arguments for matplotlib.pyplot;
        - color: string indicating color
            examples: 'b', 'g', 'r', 'c', 'm', 'y', 'k'
        - linestyle: string indicating line style
            examples: '-', '--', '-.', ':'
        - marker: string indicating marker style
            examples: 'o', '^', 's', 'D', 'v', '>', '<', 'x', '+'
        - markevery: integer indicating interval to skip marking points
            recommended: 3
            
        Returns
        -------
        None
        
        """

        if(property in shortcuts.keys()):
            property = shortcuts[property]
        # INITIATE PLOT
        if(cls.streams==[]): # new plot
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.set_xlabel('Time (y)', fontsize=16)
            ax.set_ylabel(cls.get_ylabel(property), fontsize=16)
            ax.tick_params(labelsize=14)
            cls.property = property
        else: # already plotted a stream before, axes exist
            ax = plt.gca()
            assert cls.property == property
        # MAKE AND TRACK STREAM LABEL
        cls.streams.append(strm)
        if('label' in kwargs.keys()):
            cls.labels.append(kwargs['label'])
        else:
            cls.labels.append(strm.evaluationgroup+'-'+
                              str(strm.stagenumber))
            kwargs['label'] = cls.labels[-1]+' total'
        # PLOT REQUESTED DATA
        data = getattr(strm, property)
        ax.plot(data.columns, data.sum(), **kwargs)
    
        
    @staticmethod
    def get_ylabel(property):
        """Set the y-label with units for the property plotted"""
        
        ylab = ylabel[property]
        u = unit[property]
        return(ylab.title()+' ('+u+')')
    

    @classmethod
    def reset(cls):
        """Reset the class internals to make another fresh plot"""
        
        cls.streams = []
        cls.labels = []
        cls.lines = []
        cls.color_counter = 0
        cls.color_dict = {}
        cls.marker_counter = 0
    
    
    #--------------------------------------------------------------------------
    # ELEMENTS
    #--------------------------------------------------------------------------
    @classmethod
    def add_element(cls, strm, element, **kwargs):
        """Use this function to add lines of the requested property for a
        specific element, which can be differentiated by specified color,
        marker, and line styles. Note that this may lead to changes in axes 
        limits that must be manually changed
        
        Parameters
        ----------
        cls: TimePlot class
            Class containing information about streams, plot internals, etc
        
        strm: Stream instance
            Stream or WasteForm instance containing time-dependent data for
            the property to be plotted
            
        element: str
            String indicating an element symbol abbreviation (i.e. 'cs' or 'np')
        
        kwargs: Keyword formatting arguments for matplotlib.pyplot;
        - label: string indicating data label in legend
        - color: string indicating color
            examples: 'b', 'g', 'r', 'c', 'm', 'y', 'k'
        - linestyle: string indicating line style
            examples: '-', '--', '-.', ':'
        - marker: string indicating marker style
            examples: 'o', '^', 's', 'D', 'v', '>', '<', 'x', '+'
        - markevery: integer indicating interval to skip marking points
            recommended: 3
        
        Returns
        -------
        None
        
        """
    
        ax = plt.gca()
        assert strm in cls.streams
        i = cls.streams.index(strm)
        data = getattr(cls.streams[i], cls.property)
        nuc = [n for n in data.index if
               nuclidedata.determine_nuclide_info(n)[0]==element.lower()]
        el_to_plot = sum([data.loc[n] for n in nuc])
        if('label' not in kwargs.keys()):
            kwargs['label'] = cls.labels[i]+' '+element
        ax.plot(data.columns, el_to_plot, **kwargs)

    
    #--------------------------------------------------------------------------
    # ISOTOPES
    #--------------------------------------------------------------------------
    @classmethod
    def add_isotope(cls, strm, isotope, **kwargs):
        """Add a plot for an individual isotope
        
        Parameters
        ----------
        cls: TimePlot class
            Class containing information about streams, plot internals, etc
        
        isotope: str
            String indicating isotope to add to plot (e.g. 'pu241' or 'cs137')
            
        kwargs: Keyword formatting arguments for matplotlib.pyplot;
        - label: string indicating data label in legend
        - color: string indicating color
            examples: 'b', 'g', 'r', 'c', 'm', 'y', 'k'
        - linestyle: string indicating line style
            examples: '-', '--', '-.', ':'
        - marker: string indicating marker style
            examples: 'o', '^', 's', 'D', 'v', '>', '<', 'x', '+'
        - markevery: integer indicating interval to skip marking points
            recommended: 3
            
        Returns
        -------
        None
        
        """
    
        ax = plt.gca()
        assert strm in cls.streams
        i = cls.streams.index(strm)
        data = getattr(cls.streams[i], cls.property)
        if('label' not in kwargs.keys()):
            kwargs['label'] = cls.labels[i]+' '+isotope
        ax.plot(data.columns, data.loc[isotope], **kwargs)


    #--------------------------------------------------------------------------
    # SPECIES
    # Groups of elements, such as fission products, actinides, minor actinides
    # and transuranics, can be plotted together using the species method,
    # which relies on methods in nuclidedata to match elements with groups.
    #--------------------------------------------------------------------------
    @classmethod
    def add_species(cls, strm, species, **kwargs):
        """Plot property for a single species of elements in the stream(s). 
        Probably best not to use this method if already plotting specific 
        elements or isotopes in order to avoid confusion.
        
        Parameters
        ----------
        cls: TimePlot class
        
        species: str
            Species group to add to plot; for example:
            - 'fission product' or 'fp'
            - 'actinide' or 'act'
            - 'minor actinide' or 'ma'
            - 'transuranic' or 'tru'
            - 'u3' (uraniam + protactinium)
        
        kwargs: Keyword formatting arguments for matplotlib.pyplot;
        - label: string indicating data label in legend
        - color: string indicating color
            examples: 'b', 'g', 'r', 'c', 'm', 'y', 'k'
        - linestyle: string indicating line style
            examples: '-', '--', '-.', ':'
        - marker: string indicating marker style
            examples: 'o', '^', 's', 'D', 'v', '>', '<', 'x', '+'
        - markevery: integer indicating interval to skip marking points
            recommended: 3
        
        Returns
        -------
        None
    
        """

        ax = plt.gca()
        assert strm in cls.streams
        i = cls.streams.index(strm)
        data = getattr(cls.streams[i], cls.property)
        nuc = [n for n in data.index if nuclidedata.is_group(n, species)]
        sp_to_plot = sum([data.loc[n] for n in nuc])
        if('label' not in kwargs.keys()):
            kwargs['label'] = cls.labels[i]+' '+species
        ax.plot(data.columns, sp_to_plot, **kwargs)


    #--------------------------------------------------------------------------
    # PLOT AND AXES MODIFICATIONS
    # With these methods, the user can modify axes limits, scales, labels,
    # add legends, and save and show the plots. The user can also elect to
    # hide the curves corresponding to total value for the stream property
    # that were initiated in the plot method.
    #--------------------------------------------------------------------------
    @classmethod
    def set_log(cls, *args):
        """Set axes to log scale

        Parameters
        ----------
        *args:
            'x': log scale on x-axis
            'y': log scale on y-axis

        """

        ax = plt.gca()
        if('x' in args):
            ax.set_xscale('log')
        if('y' in args):
            ax.set_yscale('log')
            ax.set_ylim(ymin=0.001)


    @classmethod
    def legend(cls, **kwargs):
        """Add a legend to the current axes object; default location is upper right
        
        Paramters
        ---------
        loc (optional): integer
            Code corresponding to location:
            0: best
            1: upper right
            2: upper left
            3: lower left
            4: lower right
            5: right
            6: center left
            7: center right
            8: lower center
            9: upper center
            10: center
            
        """
        
        ax = plt.gca()
        # just use the dang thang
        ax.legend(fontsize=16, **kwargs)


    @classmethod
    def xlim(cls, *args, **kwargs):
        """Set x-axis limits
        
        Parameters
        ----------
        *args: 
            tuple: xlim((xmin, xmax))
            values: xlim(xmin, xmax)
        
        **kwargs: Used if *args not used
            xmin: set minimum x value
            xmax: set maximum x value
        
        """

        ax = plt.gca()
        # plt.xlim should just accept these args and kwargs as is
        ax.set_xlim(*args, **kwargs)


    @classmethod
    def ylim(cls, *args, **kwargs):
        """Set y-axis limits
            
            Parameters
            ----------
            *args:
                tuple, i.e. ylim((ymin, ymax))
                values, i.e. ylim(ymin, ymax)
            
            **kwargs: Used if *args not used
                ymin: set minimum y value without changing the maximum
                ymax: set maximum y value without changing the minimum
            
            """
        
        ax = plt.gca()
        # plt.ylim should just accept these args and kwargs as is
        ax.set_ylim(*args, **kwargs)
    
    
    @classmethod
    def show(cls):
        """Show the plot"""
        
        # should add some sort of option to save first
        plt.show(block=True)
        plt.close()
        cls.reset()


    @classmethod
    def total(cls, plot=True, rescale_axes=True):
        """Do not plot the total property value, just the subvalues"""
        
        ax = plt.gca()
        if(plot==False):
            to_remove = len(cls.streams) # n 'total' lines = n streams
            ax.lines = ax.lines[to_remove:]
        if(rescale_axes):
            ax.autoscale() # for some reason, this doesn't adjust ymax
            ymax = 0.0
            for l in ax.lines:
                temp = max(l._y)
                if(temp > ymax):
                    ymax = temp
            # round to nearest order of magnitude
            ax.set_ylim(ymax = 10**(np.floor(np.log10(ymax))+1))
            

    @classmethod
    def save(cls, path=None):
        """Save the file at some place"""
        
        if(path==None): # differentiate with date and time
            path = str(datetime.date.today())+'_' # date
            path += str(datetime.datetime.now().time())[0:5]+'_'
            path += '_'+cls.property+'.png'
        plt.savefig(path)


###############################################################################


class PiePlot(object):


    """
    This class builds upon the matplotlib.pyplot functionality for fuel
    cycle plotting. In particular, it enables the plotting a snapshot of
    stream data at a specific point in time, which can be done before or
    after decay calculations have been made using the stage.decay() method. 
    
    For any stream, these methods can make pie plots using the data 
    dictionaries for isotope concentrations and heat. The decay calculations 
    load Pandas DataFrames into the stream object(s). After these have been
    loaded, the methods can utilize the Pandas Dataframes given a specific
    time point to make pie plots.
    
    Because the methods in this class alter the matplotlib.pylot internal
    state, it can be used in conjunction with matplotlib.pylot. If additional
    methods beyond those provided in the PiePlot class are required, a user
    can also import matplotlib.pyplot (referred to as 'plt' from here on)
    and, after initiating a plot using TimePlot.plot, make calls to the
    methods in plt (i.e. plt.xlim, plt.show(), plt.save(), etc).
    
    Some notes on DataFrames (fundies, if you will):
    access rows: df.loc['row name']
    access cells: df.loc['row name']['col name']
    access cols: df.loc[:]['col name']
    
    """

    # formatting defaults
    time = 0
    stream = None
    data = None
    wedges = None
    texts = None
    labels = None
    values = None
    units = ''

    #--------------------------------------------------------------------------
    # MAIN PLOT FUNCTION
    #--------------------------------------------------------------------------
    @classmethod
    def plot(cls, strm, property, time=5.0, maxslices=7, cutoff=None):
        """Generate pie chart comparing values of a Stream instance property
        
        Parameters
        ----------
        cls: PiePlot class
            Class containing information about stream, plot internals, etc
        
        strm: Stream instance
            Stream or WasteForm instance containing data for plotted property
        
        property: str
            String identical to a data attribute of the stream
            Options: 'comp', 'heat', 'radioactivity', 'radiotoxicity'; keyword
            'time' is required for 'radioactivity' and 'radiotoxicity'

        kwargs: Keyword arguments controlling plot composition
        - time: float
            Time value for data evaluation; if no value is given, method will
            plot from the 'comp' or 'heat' dictionary attributes rather than 
            Pandas DataFrame objects.
        - fraction: Boolean
            Indicates whether to plot fractions rather than absolute values
        - cutoff: float
            Minimum below which values are not included in plot
        - maxslices: int
            Maximum number of species to include in plot

        Returns
        -------
        None; generates matplotlib figure
        
        """
        
        if(property in shortcuts.keys()):
            property = shortcuts[property]
        data = cls._get_data_dict(strm, property, time)
        labels, values = cls._sort_data(data)
        # apply maxslices; flip order; leave room for "other" slice
        labels, values = cls._apply_maxslices(labels, values, maxslices)
        # apply cutoff, make fraction, add slice for "other"
        labels, values = cls._apply_cutoff(labels, values, data, cutoff)
        # plot what ya got
        fig, ax = plt.subplots(figsize=(10, 8),subplot_kw=dict(aspect="equal"))
        # assign to class variables for additional methods
        cls.wedges, cls.texts = ax.pie(values, startangle=-40)
        #wedgeprops=dict(width=0.5)
        cls.labels, cls.values = [labels, values]
        cls.units = cls.get_units(strm, property)
        cls.data = data
    
    
    @staticmethod
    def get_units(strm, property):
        """Set the y-label with units for the property plotted"""
        
        try: # get units based on property
            units = strm.units.keys()[strm.units.values().index(property)]
        except AttributeError:
            units = ''
        return(units)
    
    
    #--------------------------------------------------------------------------
    # DATA RETRIEVAL FROM DF/DICT
    #--------------------------------------------------------------------------
    @classmethod
    def _get_data_dict(cls, strm, property, time):
        """Retrieve dictionary from stream containing data to plot"""
        
        data = cls._interpolate_df(strm, property, time)
        data = data.to_dict()
        return(data)


    @staticmethod
    def _interpolate_df(strm, property, target):
        """Get values for all nuclides for a specific OPUS output at a
        specified time using interpolation over pandas df to plot pie chart
        
        Parameters
        ----------
        strm: Stream instance
            Must contain Pandas DF corresponding to requested property
        
        property: str
            comp, heat, radioactivity, radiotoxicity
        
        target: float
            Time to return interpolated values
        
        Returns
        -------
        Pandas Series object at given time with interpolated values
        for all nuclides
        
        """
        
        temp = getattr(strm, property)
        if(target in temp.columns):
            return(temp[target])
        t1 = min(temp.columns, key=lambda x: abs(float(x)-target))
        idx = temp.columns.get_loc(t1)
        t2 = temp.columns[int(-np.sign(float(t1)-target)+idx)]
        temp = temp.loc[:, [t1, t2]]
        new_cols = {}
        for x in temp.columns:
            new_cols[x]=float(x)
        temp = temp.rename(columns=new_cols)
        temp.insert(1, target, pd.Series(np.nan), allow_duplicates=True)
        temp = temp.sort_index(axis=1)
        temp = temp.astype(float)
        temp = temp.apply(lambda x: x.interpolate(method='index'), axis=1)
        return(temp[temp.columns[1]])


    #--------------------------------------------------------------------------
    # DATA MANIPULATION AND SORTING
    #--------------------------------------------------------------------------
    @staticmethod
    def _sort_data(data_dict):
        """Given a dictionary, sort labels, values in descending order"""

        l = data_dict.keys()
        v = data_dict.values()
        l = [x for _,x in sorted(zip(v,l))]
        v = sorted(v)
        return(l, v)
    
    
    @staticmethod
    def _apply_maxslices(labels, values, maxslices):
        """Apply maxslices cutoff"""
        
        v = values[len(values)-maxslices:][::-1]
        l = labels[len(labels)-maxslices:][::-1]
        return(l, v)
    
    
    @staticmethod
    def _apply_cutoff(labels, values, data_dict, cutoff):
        """"Apply cutoff, make a slice in the pie to account for 
        other values not shown"""
        
        v = [val for val in values if val > cutoff]
        l = labels[len(labels)-len(v):]
        # add slice for "other"
        l = l+['Other']
        v = v+[sum(data_dict.values())-sum(v)]
        return(l, v)


    #--------------------------------------------------------------------------
    # LABELS AND ANNOTATION
    #--------------------------------------------------------------------------
    @classmethod
    def label(cls):
        """Add slice labels to pie plot"""
        
        ax = plt.gca()
        plot_labels = cls._make_plot_labels()
        bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
        kw = dict(xycoords='data', textcoords='data',
                  arrowprops=dict(arrowstyle="-"),
                  bbox=bbox_props, zorder=1, va="center")
        for i, p in enumerate(cls.wedges):
            ang = (p.theta2 - p.theta1)/2. + p.theta1
            y = np.sin(np.deg2rad(ang))
            x = np.cos(np.deg2rad(ang))
            horizontalalignment = "center"
            #{-1: "right", 1: "left"}[int(np.sign(x))]
            connectionstyle = "angle,angleA=0,angleB={}".format(ang)
            kw["arrowprops"].update({"connectionstyle": connectionstyle})
            ax.annotate(plot_labels[i], xy=(x, y), xytext=(1.4*x, 1.4*y),
                        horizontalalignment=horizontalalignment, **kw)


    @classmethod
    def _make_plot_labels(cls):
        """Make some annotations for the figure"""
        
        plot_labels = []
        for i in range(0, len(cls.labels)):
            l = cls.labels[i].title()+'\n'
            v = '%.2f'%Decimal(cls.values[i]*100.0/sum(cls.data.values()))
            l += v+'%'+'\n'
            l += str(round(cls.values[i], 2))+' '+cls.units
            plot_labels.append(l)
        return(plot_labels)


    @classmethod
    def annotate(cls):
        """Add annotations for total heat and time evaluated"""
        return()
    
    
    #--------------------------------------------------------------------------
    # DISPLAY
    #--------------------------------------------------------------------------
    @classmethod
    def show(cls):
        """Show the plot"""
        
        # should add some sort of option to save first
        plt.show(block=True)
        plt.close()
        cls.reset()


    @classmethod
    def reset(cls):
        """Reset the class internals to make another fresh plot"""
        
        cls.streams = []
        cls.labels = []
        cls.lines = []
        cls.color_counter = 0
        cls.color_dict = {}
        cls.marker_counter = 0


# TO ADD
# annotate
# legend
# modify "label" to take formatting arguments


###############################################################################


# lots of colors
# https://stackoverflow.com/questions/8389636/creating-over-20-unique-legend-colors-using-matplotlib
#
# labeling pies
# https://matplotlib.org/gallery/pie_and_polar_charts/pie_and_donut_labels.html



###############################################################################


class LoadPlot(object):


    """
    The methods in this object plot the inequality and equality 
    constraints that bound the solution space for waste loading. 
    If possible, the plot of that solution space can be generated 
    when loading a waste form.

    """
    
    @classmethod
    def plot(cls, a_ub, b_ub, key_ub, a_eq=None, b_eq=None,
             key_eq=None, x_max=None, label=True, units=''):
        """Plot the waste loading solution space

        Parameters
        ----------
        cls: LoadPlot class object
        
        a_ub: (M,2) array-like
            A matrix containing upper bound inequality constraints

        b_ub: (M,1) array-like
            Right-hand side matrix in Ax <= b
        
        key_ub: list
            List of strings describing the upper bound constraints

        a_eq (optional): (M,2) array-like
            A matrix containing upper bound inequality constraints

        b_eq (optional): numpy array
            Right-hand side matrix in Ax = b
            
        key_eq (optional): list
            List of strings describing the equality constraints

        x_max (optional): float
            Maximum value for the x-axis

        label (optional): Boolean
            Indicates whether to label the constraint lines on the figure

        Returns
        -------
        None, but displays a matplotlib figure

        """
        
        a_ub, b_ub, key_ub = cls._preprocess_ineq(a_ub, b_ub, key_ub)
        x_plt = cls._xbound(b_ub, x_max)
        fig, ax = plt.subplots()
        # generate ineq lines
        cls._generate_lines(x_plt, a_ub, b_ub, key_ub)
        # generate eq lines
        if(not any([inp==None for inp in [a_eq, b_eq, key_eq]])):
            cls._generate_lines(x_plt, a_ub, b_ub, key_ub)
        if(label):
            LabelLines.label_lines(plt.gca().get_lines(), zorder=2.5)
        ax.set_xlim([0, max(x_plt)])
        ax.set_ylim([0, max(x_plt)])
        if(units != ''):
            units = ' ('+units+')'
        ax.set_xlabel('Mass of diluent matrix'+units)
        ax.set_ylabel('Mass of waste'+units)
        plt.show(block=True)
        plt.close()


    @staticmethod
    def _preprocess_ineq(a, b, key):
        """ """
        
        to_delete = []
        # skip the nonzero constraints
        for i in range(0, len(key)):
            if('nonzero' in key[i].lower()):
                to_delete.append(i)
        a = np.delete(a, to_delete, axis=0)
        b = np.delete(b, to_delete)
        key = [k for k in key if 'nonzero' not in k.lower()]
        return(a, b, key)
        
    
    @staticmethod
    def _xbound(b, xmax):
        """ """
        
        if(xmax):
            xplt = np.arange(0.0, xmax, 20.0)
            xplt = np.append(xplt, xmax)
        else:
            xplt = np.arange(0.0, b[0], 20.0)
            xplt = np.append(xplt, b[0])
        return(xplt)


    @classmethod
    def _generate_lines(cls, xpts, a, b, key, ineq=True):
        """ """
        
        y_lb = np.zeros(len(xpts)) # set lower bound
        y_ub = None # set by values
        ax = plt.gca()
        for row in range(0, len(b)):
            ypts = np.zeros(len(xpts))
            for x in range(0, len(xpts)):
                ypts[x] = (b[row]-a[row][1]*xpts[x])/a[row][0]
            #temp = temp/1e3 # /1e3: g->kg
            ax.plot(xpts, ypts, label=key[row])
            if(ineq):
                y_lb,y_ub = cls._shade_ineq(xpts, ypts, a[row], y_lb, y_ub)
        if(ineq):
            ax.fill_between(xpts, y_lb, y_ub, where=y_lb<y_ub,
                            facecolor='gray', alpha=0.5)


    @staticmethod
    def _shade_ineq(x, y, c, ylb, yub):
        """Update shading of solution region"""
        
        if(c[0]<0.0): # "greater than"
            ylb = np.maximum(ylb, y)
        else: # "less than"
            if(yub is not None):
                yub = np.minimum(yub, y)
            else: # upper bound not defined yet
                yub = y
        return(ylb, yub)


###############################################################################


class LabelLines(object):


    """
    Line-labeling module developed by StackOverflow user NauticalMile; 
    modified into Python class by Milos Atz.
    
    https://stackoverflow.com/questions/16992038/inline-labels-in-matplotlib
    
    NOTE (2018-04-25): 
    Want to ignore depreciation warning associated with ax = line.get_axes()
    (see above, in import statements)
    
    """


    @classmethod
    def label_lines(cls, lines, align=True, xvals=None, **kwargs):
        """Label all lines in a figure"""
        
        ax = lines[0].get_axes()
        lab_lines = []
        labels = []
        for line in lines: # only lines w/ labels other than default
            label = line.get_label()
            if "_line" not in label:
                lab_lines.append(line)
                labels.append(label)
        if xvals is None:
            xmin,xmax = ax.get_xlim()
            xvals = np.linspace(xmin,xmax,len(lab_lines)+2)[1:-1]
        for line,x,label in zip(lab_lines,xvals,labels):
            cls.label_line(line,x,label,align,**kwargs)


    @staticmethod
    def label_line(line, x, label=None, align=True, **kwargs):
        """Label line with line2D label data"""

        ax = line.get_axes()
        xdata = line.get_xdata()
        ydata = line.get_ydata()
        if (x < xdata[0]) or (x > xdata[-1]):
            print('x label location is outside data range!')
            return
        # Find corresponding y co-ordinate and angle of the line
        ip = 1
        for i in range(len(xdata)):
            if x < xdata[i]:
                ip = i
                break
        y = ydata[ip-1] + (ydata[ip]-ydata[ip-1])*(x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])
        if(y > max(xdata)):
            ip = 1
            y = max(xdata)-800.0
            for i in range(len(xdata)):
                if(y < ydata[i]):
                    ip = i
                    break
        x = xdata[ip]
        if not label:
            label = line.get_label()
        if align: # Compute the slope
            dx = xdata[ip] - xdata[ip-1]
            dy = ydata[ip] - ydata[ip-1]
            ang = degrees(atan2(dy,dx))
            if(ang<0):
                trans_angle = ang+4.0
            else:
                trans_angle = ang-4.0
        else:
            trans_angle = 0
        #Set a bunch of keyword arguments
        if 'color' not in kwargs:
            kwargs['color'] = line.get_color()
        if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
            kwargs['ha'] = 'center'
        if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
            kwargs['va'] = 'center'
        if 'backgroundcolor' not in kwargs:
            kwargs['backgroundcolor'] = ax.get_axis_bgcolor()
        if 'clip_on' not in kwargs:
            kwargs['clip_on'] = True
        if 'zorder' not in kwargs:
            kwargs['zorder'] = 2.5
        ax.text(x,y,label,rotation=trans_angle,**kwargs)


###############################################################################


