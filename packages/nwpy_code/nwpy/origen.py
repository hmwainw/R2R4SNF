###############################################################################
"""
    Last edited on July 13, 2019
    
    @author: matz
    
    comments: Controls the manipulation of data to make ORIGEN-S input files.
    
"""
###############################################################################
import os
import copy
import datetime
import subprocess
import warnings
import pandas as pd
import copy
from decimal import Decimal
from nwpy import nuclidedata
mass_cutoff = 1e-9 # cutoff for excluding nuclides from ORIGEN input
###############################################################################


class Origen:
    
    
    """
    The main purpose of this code is to write, run, and read ORIGEN decay 
    cases for the material streams in the FCES analysis examples. In this 
    section, the code that controls the input file generation and call to 
    ORIGEN are written. In additional files (build.py, template.py) are 
    the general functions that are used to write the ORIGEN files.

    """
    
    
    def __init__():
        pass
    
    
    @classmethod
    def make_file(cls, wdpath, name, strm_inst, starttime, endtime, steps,
                  interpolation, cooling=False):
        """Collect the information necessary to build the ORIGEN input
        file from the relevant stream and build it using the methods 
        in the this module.
        
        Parameters
        ----------
        wdpath: str
            Path to working directory for file I/O
        
        name: str
            Name of Stage instance
            
        strm_inst: Stream instance
            Stream for which the ORIGEN input file should be generated
        
        starttime: float
            Time at which composition data is read; must reflect value
            in column names of .comp DataFrame attribute.
            
        endtime: float
            End time for ORIGEN decay calculation
        
        steps: int
            Number of steps to be used for the ORIGEN calculation
        
        interpolation: str
            Indicates step interpolation method
            'I': Linear
            'L': Logarithmic
        
        cooling: boolean, optional
            If true, input will simulate 5 y cooling before recycle.
        
        Returns
        -------
        Filename of newly written input
        
        """
        
        # generate and parse the time info
        time_input = format_time(endtime, steps, interpolation)
        # format the nuclide material info
        temp = strm_inst.comp[starttime].to_dict() # cast as dict
        nuclide_input = format_materials(temp)
        case = cls.make_case(strm_inst, name)
        if(cooling):
            if(endtime==0.0):
                saved = '[0]'
            else:
                saved = 'LAST'
            case += '_cool'
            output = ['grams', 'watts']#, 'curies', 'H2OM**3']
        else:
            saved = 'ALL' # time steps where data is saved
            output = ['grams', 'watts', 'curies', 'H2OM**3'] # opus outputs
        path = case+'.f71'
        file = case+'.inp'
        t = write_temp(nuclide_input, time_input, path, saved, output)
        a = open(os.path.join(wdpath, file), 'w')
        a.write(t)
        a.close()
        return(file)
    
    
    @staticmethod
    def make_case(strm, stg_name):
        """Make the file name for the ORIGEN calculation based on attributes 
        of the Stream or WasteForm instance used as input"""
    
        #name = strm.evaluationgroup+'-'+str(strm.stagenumber) # egXX-X
        name = stg_name+'_'+strm.form #
        if(hasattr(strm, 'batch') and 'batch' not in strm.form):
            name += '_batch'+str(strm.batch)
        if(hasattr(strm, 'loading') and strm.id=='snf'):
            name += '_'+str(strm.loading)
        if(hasattr(strm, 'loading_fraction') and strm.id=='hlw'):
            name += '_'+str(strm.loading_fraction)
        return(name)
    
    
    @classmethod
    def run_origen(cls, wdpath, infile):
        """Run an ORIGEN input file
        
        Parameters
        ----------
        self: Stage instance
        
        wdpath: str
            Path to working directory where input file is located
            
        infile: str
            File name ('.inp') of ORIGEN input file
        
        Returns
        -------
        None
        
        """
        
        # Note: SCALE must be added to PATH! Can confirm via test
        path_to_scale = cls._find_scale()
        assert os.path.exists(path_to_scale),'Path to scalerte executable not found'
        path_to_infile = os.path.join(wdpath, infile)
        subprocess.call([path_to_scale, path_to_infile])


    @staticmethod
    def _find_scale():
        """Find the path to the scalerte executable based on the contents 
        of the .bashrc path environment variables"""
        
        path_to_scale = []
        for p in os.environ['PATH'].split(os.pathsep):
            try: # if not a directory, this will fail
                temp = os.listdir(p)
            except OSError:
                continue
            else:
                if('scalerte' in os.listdir(p)): # unix
                    path_to_scale.append(os.path.join(p, 'scalerte'))
                elif('scalerte.exe' in os.listdir(p)): # windows
                    path_to_scale.append(os.path.join(p, 'scalerte.exe'))
                else:
                    continue
        warning = ('Path to scalerte executable not found. Make sure '+
                   'it is added to .bashrc PATH environment')
        assert path_to_scale != [], warning
        assert os.path.exists(path_to_scale[0]), warning
        return(path_to_scale[0])
        
    
    @staticmethod
    def update_stream(wdpath, str_inst, infile, starttime):
        """Update Stream attributes (Pandas DataFrames) that contain data
        for composition, heat radioactivity and radiotoxicity from OPUS 
        plt files.
        
        Parameters
        ----------
        self: Stage instance
        
        str_inst: Stream instance
            Stream to be updated
        
        infile: str
            File name of the ORIGEN calculation
        
        starttime: float
        
        Returns
        -------
        New Stream instance with new and/or updated attributes
        
        """
        
        strm = copy.deepcopy(str_inst)
        # unit/attribute associations
        u = strm.units
        # get plt files first
        pltfiles = [os.path.join(wdpath, x) for x in os.listdir(wdpath)
                    if 'plt' in x and infile[0:len(infile)-3] in x]
        for file in pltfiles: # read opus data
            txtlines = open(file).read().splitlines()
            unit = txtlines[2]
            attr = u[unit] # determine attribute based on units
            dat=pd.read_csv(file, delim_whitespace=True,header=4,index_col=0)
            # remove duplicate columns
            dat=dat.loc[:,~dat.columns.str.replace("(\.\d+)$","").duplicated()]
            dat.columns = pd.to_numeric(dat.columns)
            dat = dat.drop(['subtotal','total'])
            strm.update(attr, dat, starttime)
        return(strm)


###############################################################################
# MODULE LEVEL FUNCTIONS
###############################################################################


def format_materials(comp_dict):
    """Create a list of isotopes and masses for ORIGEN input

    Parameters
    ----------
    comp_dict: dict
        Dictionary with keys: nuclides and values: masses
    
    Returns
    -------
    str with the materials and values formatted for ORIGEN input
    
    """
    
    comp = copy.deepcopy(comp_dict)
    comp = nuclidedata.group_daughters(pd.Series(comp)).to_dict()
    mat_list = []
    for key in comp:
        if(comp[key] > mass_cutoff):
            mat_list.append(key.strip()+'='+'%.5e' % Decimal(comp[key]))
    return('[ '+' '.join(mat_list)+' ]')


def format_time(endtime, steps, interpolation):
    """Create the decay time input information for ORIGEN
    
    Parameters
    ----------
    endtime: float
        End time for ORIGEN decay calculation
    
    steps: int
        Number of steps to be used for the ORIGEN calculation
    
    interpolation: str
        Indicates step interpolation method; 'I': Linear, 'L': Logarithmic
    
    Returns
    -------
    str with formatted time input for ORIGEN-S input file
    
    """
    
    if(endtime==0.0): # on savio, timestep must be positive; use a small
        return('0.1') # timestep, then have opus return starting values
    else:
        if(steps > 0.0):
            steps = int(steps-1)
        else:
            raise ValueError('Number of steps must be greater than 0')
            #steps = int(steps) # i.e. n=0
        return('['+str(steps)+interpolation+' 1 '+str(endtime)+']')
        

def format_opus_nuclides(nuclide_str):
    """Organize the list of nuclides that are to be included in the OPUS input
    
    Parameters
    ----------
    nuclide_str: str
        Space-separated nuclides to be included in the OPUS input
        
    Returns
    -------
    1. Nuclide list formatted for writing OPUS input
    2. Number of nuclides included in OPUS input
    
    """
    nuclide_list = nuclide_str.split(' ')
    formatted_nuclide_str = ''
    sublist = []
    for nuc in nuclide_list:
        sublistlen = sum([len(x) for x in sublist])+len(sublist)
        if(sublistlen+len(nuc)+1 < 72):
            sublist.append(nuc)
        else:
            formatted_nuclide_str += ' '+' '.join(sublist)+'\n'
            sublist = []
            sublist.append(nuc)
    # catch the last one; there will always be a last one!
    formatted_nuclide_str += ' '+' '.join(sublist)
    return(formatted_nuclide_str, len(nuclide_list))


def write_temp(matlist, timelist, datapath, savedtimes, units):
    """Assemble the ORIGEN input file from templates and replace keywords
    
    Parameters
    ----------
    matlist: str
        List of nuclides and their masses
    
    timelist: str
        Time information for ORIGEN calculation
    
    datapath: str
        Path to the ORIGEN binary output file for OPUS plotting
    
    savedtimes: str
        Either 'LAST' or 'ALL' depending on the type of calculation
    
    units: str
        Units for data to be plotted by OPUS
    
    Returns
    -------
    str containing written ORIGEN input file
    
    """
    temp = preamble_temp
    temp = temp.replace("TAG-DATETIME", "%s" %(datetime.datetime.now()))
    temp = temp+decay_case_temp
    temp = temp.replace("TAG-NUCLIDEINFO", "%s" %(matlist))
    temp = temp.replace("TAG-TIMEINFO", "%s" %(timelist))
    temp = temp.replace("TAG-DATAPATH", "%s" %(datapath))
    temp = temp.replace("TAG-SAVETIMESTEPS", "%s" %(savedtimes))
    # get the nuclides to plot
    nuclides_to_plot, nrank = format_opus_nuclides(nuclides)
    for plt in range(0,len(units)):
        temp = temp+opus_temp
        temp = temp.replace("TAG-DATAPATH", "%s" %(datapath))
        temp = temp.replace("TAG-OUTPUTUNITS", "%s" %(units[plt]))
        temp = temp.replace("TAG-NRANK", "%s" %(nrank))
        temp = temp.replace("TAG-SYMNUC", "%s" %(nuclides_to_plot))
    return(temp)


###############################################################################
# TEMPLATES AND DATA FOR ORIGEN-S INPUTS
###############################################################################
# PREAMBLE TEMPLATE
#------------------------------------------------------------------------------
preamble_temp="""'WRITTEN BY SCALE6.2.2 ORIGEN-S INPUT GENERATOR
'PART OF THE NWPY PYTHON PACKAGE
'MILOS ATZ
'WASTE MANAGEMENT GROUP
'DEPARTMENT OF NUCLEAR ENGINEERING
'UNIVERSITY OF CALIFORNIA BERKELEY
'Code executed: TAG-DATETIME
'
=origen
"""
#------------------------------------------------------------------------------
# DECAY CASE TEMPLATECASE BY CASE TEMPLATE
#------------------------------------------------------------------------------
decay_case_temp="""
case(decay){
    
    lib{ file="end7dec" }
    
    mat{
        units=GRAMS
        iso=TAG-NUCLIDEINFO
    }
    
    time{
        units=YEARS
        t=TAG-TIMEINFO
    }
    save{ file="TAG-DATAPATH" steps=TAG-SAVETIMESTEPS }
    }
end
"""
#------------------------------------------------------------------------------
# OPUS PLOT TEMPLATE
#------------------------------------------------------------------------------
# This is the plotting template that will be inserted at the end of the input
# file. It's much easier to parse the output data when ORIGEN plots it for
# you. In the output, the data is arranged into neat data tables that will
# be easy to import into python. Additionally, plotting in OPUS is the
# straightforward way to get non-mass quantities like heat, radioactivity,
# and radiotoxicity.
#
# Note: The OPUS plotting is not split by case. So I can plot all of the
# saved data points in a single plot command - this should not cause the
# the command to overrun the ORIGEN character limit until about 20 time
# steps; 10 billion years is 18 time steps.
#------------------------------------------------------------------------------
opus_temp="""
=opus
data="TAG-DATAPATH"
units=TAG-OUTPUTUNITS
time=year
libtype=all
nrank=TAG-NRANK
symnuc =TAG-SYMNUC end
end
"""
#------------------------------------------------------------------------------
# NUCLIDE LIST
#------------------------------------------------------------------------------
# Nuclides to be included in the OPUS plotter; this represents *all* the
# nuclides in the ORIGEN: actinide, fission product, and activation product
# libraries. That means that many will be zero.
#------------------------------------------------------------------------------
nuclides = ''
# hydrogen
nuclides += 'H-1 H-2 H-3 '
# helium
nuclides += 'HE-3 HE-4 HE-5 HE-6 '
# lithium
nuclides += 'LI-6 LI-7 LI-8 '
# beryllium
nuclides += 'BE-7 BE-8 BE-9 BE-10 BE-11 '
# boron
nuclides += 'B-10 B-11 B-12 '
# carbon
nuclides += 'C-12 C-13 C-14 C-15 '
# nitrogen
nuclides += 'N-13 N-14 N-15 N-16 '
# oxygen
nuclides += 'O-16 O-17 O-18 O-19 '
# fluorine
nuclides += 'F-19 F-20 '
# neon
nuclides += 'NE-20 NE-21 NE-22 NE-23 '
# sodium
nuclides += 'NA-22 NA-23 NA-24 NA-24M NA-25 '
# magnesium
nuclides += 'MG-24 MG-25 MG-26 MG-27 MG-28 '
# aluminum
nuclides += 'AL-26 AL-27 AL-28 AL-29 AL-30 '
# silicon
nuclides += 'SI-28 SI-29 SI-30 SI-31 SI-32 '
# phosphorus
nuclides += 'P-31 P-32 P-33 P-34 '
# sulfur
nuclides += 'S-32 S-33 S-34 S-35 S-36 S-37 '
# chlorine
nuclides += 'CL-35 CL-36 CL-37 CL-38 CL-38M '
# argon
nuclides += 'AR-36 AR-37 AR-38 AR-39 AR-40 AR-41 AR-42 '
# potassium
nuclides += 'K-39 K-40 K-41 K-42 K-43 K-44 '
# calcium
nuclides += 'CA-40 CA-41 CA-42 CA-43 CA-44 CA-45 CA-46 CA-47 CA-48 CA-49 '
# scandium
nuclides += 'SC-44 SC-44M SC-45 SC-45M SC-46 SC-46M SC-47 SC-48 SC-49 SC-50 '
# titanium
nuclides += 'TI-44 TI-45 TI-46 TI-47 TI-48 TI-49 TI-50 TI-51 '
# vanadium
nuclides += 'V-48 V-49 V-50 V-51 V-52 V-53 V-54 '
# chromium
nuclides += 'CR-48 CR-49 CR-50 CR-51 CR-52 CR-53 CR-54 CR-55 CR-66 CR-67 '
# manganese
nuclides += 'MN-52 MN-53 MN-54 MN-55 MN-56 MN-57 '+\
            'MN-58 MN-66 MN-67 MN-68 MN-69 '
# iron
nuclides += 'FE-54 FE-55 FE-56 FE-57 FE-58 FE-59 FE-60 FE-65 FE-66 '+\
            'FE-67 FE-68 FE-69 FE-70 FE-71 FE-72 '
# cobalt
nuclides += 'CO-55 CO-56 CO-57 CO-58M CO-58 CO-59 CO-60 CO-60M CO-61 '+\
            'CO-62 CO-65 CO-66 CO-67 CO-68 CO-69 CO-70 CO-71 CO-72 '+\
            'CO-73 CO-74 CO-75 '
# nickel
nuclides += 'NI-56 NI-57 NI-58 NI-59 NI-60 NI-61 NI-62 NI-63 NI-64 '+\
            'NI-65 NI-66 NI-67 NI-68 NI-69 NI-70 NI-71 NI-72 NI-73 '+\
            'NI-74 NI-75 NI-76 NI-77 NI-78 '
# copper
nuclides += 'CU-62 CU-63 CU-64 CU-65 CU-66 CU-67 CU-68 CU-68M CU-69 '+\
            'CU-70 CU-70M CU-71 CU-72 CU-73 CU-74 CU-75 CU-76 CU-77 '+\
            'CU-78 CU-79 CU-80 '
# zinc
nuclides += 'ZN-63 ZN-64 ZN-65 ZN-66 ZN-67 ZN-68 ZN-69 ZN-69M ZN-70 '+\
            'ZN-71 ZN-71M ZN-72 ZN-73 ZN-74 ZN-75 ZN-76 ZN-77 ZN-78 '+\
            'ZN-79 ZN-80 ZN-81 ZN-82 ZN-83 '
# gallium
nuclides += 'GA-67 GA-68 GA-69 GA-70 GA-71 GA-72 GA-72M GA-66 GA-73 '+\
            'GA-74 GA-74M GA-75 GA-76 GA-77 GA-78 GA-79 GA-80 GA-81 '+\
            'GA-82 GA-83 GA-84 GA-85 GA-86 '
# germanium
nuclides += 'GE-66 GE-67 GE-68 GE-69 GE-70 GE-71 GE-71M GE-72 GE-73 '+\
            'GE-73M GE-74 GE-75 GE-75M GE-76 GE-77 GE-77M GE-78 GE-79 '+\
            'GE-79M GE-80 GE-81 GE-81M GE-82 GE-83 GE-84 GE-85 GE-86 '+\
            'GE-87 GE-88 GE-89 '
# arsenic
nuclides += 'AS-69 AS-71 AS-72 AS-73 AS-74 AS-75 AS-75M AS-76 AS-77 '+\
            'AS-78 AS-79 AS-80 AS-81 AS-82 AS-82M AS-83 AS-84 AS-85 '+\
            'AS-86 AS-87 AS-88 AS-89 AS-90 AS-91 AS-92 '
# selenium
nuclides += 'SE-72 SE-73 SE-73M SE-74 SE-75 SE-76 SE-77 SE-77M SE-78 '+\
            'SE-79 SE-79M SE-80 SE-81 SE-81M SE-82 SE-83 SE-83M SE-84 '+\
            'SE-85 SE-86 SE-87 SE-88 SE-89 SE-90 SE-91 SE-92 SE-93 SE-94 '
# bromine
nuclides += 'BR-75 BR-76 BR-77 BR-77M BR-78 BR-79 BR-79M BR-80 BR-80M '+\
            'BR-81 BR-82 BR-82M BR-83 BR-84 BR-84M BR-85 BR-86 BR-87 '+\
            'BR-88 BR-89 BR-90 BR-91 BR-92 BR-93 BR-94 BR-95 BR-96 BR-97 '
# krypton
nuclides += 'KR-76 KR-77 KR-78 KR-79 KR-79M KR-80 KR-81 KR-81M KR-82 '+\
            'KR-83 KR-83M KR-84 KR-85 KR-85M KR-86 KR-87 KR-88 KR-89 '+\
            'KR-90 KR-91 KR-92 KR-93 KR-94 KR-95 KR-96 KR-97 KR-98 '+\
            'KR-99 KR-100 '
# rubidium
nuclides += 'RB-79 RB-81 RB-82 RB-83 RB-84 RB-85 RB-86 RB-86M RB-87 '+\
            'RB-88 RB-89 RB-90 RB-90M RB-91 RB-92 RB-93 RB-94 RB-95 '+\
            'RB-96 RB-97 RB-98 RB-99 RB-100 RB-101 RB-102 '
# stronitum
nuclides += 'SR-82 SR-83 SR-84 SR-85 SR-85M SR-86 SR-87 SR-87M SR-88 '+\
            'SR-89 SR-90 SR-91 SR-92 SR-93 SR-94 SR-95 SR-96 SR-97 '+\
            'SR-98 SR-99 SR-100 SR-101 SR-102 SR-103 SR-104 SR-105 '
# yttrium
nuclides += 'Y-85 Y-86 Y-87 Y-87M Y-88 Y-89 Y-89M Y-90 Y-90M Y-91 '+\
            'Y-91M Y-92 Y-93 Y-93M Y-94 Y-95 Y-96 Y-96M Y-97 Y-97M '+\
            'Y-98 Y-98M Y-99 Y-100 Y-101 Y-102 Y-103 Y-104 Y-105 '+\
            'Y-106 Y-107 Y-108 '
# zirconium
nuclides += 'ZR-86 ZR-87 ZR-88 ZR-89 ZR-89M ZR-90 ZR-90M ZR-91 ZR-92 '+\
            'ZR-93 ZR-94 ZR-95 ZR-96 ZR-97 ZR-98 ZR-99 ZR-100 ZR-101 '+\
            'ZR-102 ZR-103 ZR-104 ZR-105 ZR-106 ZR-107 ZR-108 ZR-109 '+\
            'ZR-110 '
# niobium
nuclides += 'NB-89 NB-90 NB-90M NB-91 NB-91M NB-92 NB-92M NB-93 NB-93M '+\
            'NB-94 NB-94M NB-95 NB-95M NB-96 NB-97 NB-97M NB-98 NB-98M '+\
            'NB-99 NB-99M NB-100 NB-100M NB-101 NB-102 NB-102M NB-103 '+\
            'NB-104 NB-104M NB-105 NB-106 NB-107 NB-108 NB-109 NB-110 '+\
            'NB-111 NB-112 NB-113 '
# molybdenum
nuclides += 'MO-90 MO-91 MO-92 MO-93M MO-93 MO-94 MO-95 MO-96 MO-97 '+\
            'MO-98 MO-99 MO-100 MO-101 MO-102 MO-103 MO-104 MO-105 '+\
            'MO-106 MO-107 MO-108 MO-109 MO-110 MO-111 MO-112 MO-113 '+\
            'MO-114 MO-115 '
# technetium
nuclides += 'TC-93 TC-95 TC-95M TC-96 TC-97 TC-97M TC-98 TC-99 TC-99M '+\
            'TC-100 TC-101 TC-102 TC-102M TC-103 TC-104 TC-105 TC-106 '+\
            'TC-107 TC-108 TC-109 TC-110 TC-111 TC-112 TC-113 TC-114 '+\
            'TC-115 TC-116 TC-117 TC-118 '
# ruthenium
nuclides += 'RU-95 RU-96 RU-97 RU-98 RU-99 RU-100 RU-101 RU-102 RU-103 '+\
            'RU-104 RU-105 RU-106 RU-107 RU-108 RU-109 RU-110 RU-111 '+\
            'RU-112 RU-113 RU-114 RU-115 RU-116 RU-117 RU-118 RU-119 '+\
            'RU-120 '
# rhodium
nuclides += 'RH-99 RH-99M RH-100 RH-101 RH-101M RH-102 RH-102M RH-103 '+\
            'RH-103M RH-104 RH-104M RH-105 RH-105M RH-106 RH-106M '+\
            'RH-107 RH-108 RH-108M RH-109 RH-110 RH-110M RH-111 RH-112 '+\
            'RH-113 RH-114 RH-115 RH-116 RH-117 RH-118 RH-119 RH-120 '+\
            'RH-121 RH-122 '
# palladium
nuclides += 'PD-99 PD-100 PD-101 PD-102 PD-103 PD-104 PD-105 PD-106 '+\
            'PD-107 PD-107M PD-108 PD-109 PD-109M PD-110 PD-111 '+\
            'PD-111M PD-112 PD-113 PD-114 PD-115 PD-116 PD-117 '+\
            'PD-118 PD-119 PD-120 PD-121 PD-122 PD-123 PD-124 '
# silver
nuclides += 'AG-103 AG-105 AG-105M AG-106 AG-106M AG-107 AG-107M AG-108 '+\
            'AG-108M AG-109 AG-109M AG-110 AG-110M AG-111 AG-111M AG-112 '+\
            'AG-113 AG-113M AG-114 AG-115 AG-115M AG-116 AG-116M AG-117 '+\
            'AG-117M AG-118 AG-118M AG-119 AG-120 AG-120M AG-121 AG-122 '+\
            'AG-122M AG-123 AG-124 AG-125 AG-126 AG-127 AG-128 AG-129 AG-130 '
# cadmium
nuclides += 'CD-105 CD-106 CD-107 CD-108 CD-109 CD-110 CD-111 CD-111M '+\
            'CD-112 CD-113 CD-113M CD-114 CD-115 CD-115M CD-116 CD-117 '+\
            'CD-117M CD-118 CD-119 CD-119M CD-120 CD-121 CD-121M CD-122 '+\
            'CD-123 CD-123M CD-124 CD-125 CD-126 CD-127 CD-128 CD-129 '+\
            'CD-130 CD-131 CD-132 '
# indium
nuclides += 'IN-107 IN-109 IN-111 IN-111M IN-112 IN-112M IN-113 IN-113M '+\
            'IN-114 IN-114M IN-115 IN-115M IN-116 IN-116M IN-117 IN-117M '+\
            'IN-118 IN-118M IN-119 IN-119M IN-120 IN-120M IN-121 IN-121M '+\
            'IN-122 IN-122M IN-123 IN-123M IN-124 IN-124M IN-125 IN-125M '+\
            'IN-126 IN-126M IN-127 IN-127M IN-128 IN-128M IN-129 IN-129M '+\
            'IN-130 IN-130M IN-131 IN-131M IN-132 IN-133 IN-134 IN-135 '
# tin
nuclides += 'SN-111 SN-112 SN-113 SN-113M SN-114 SN-115 SN-116 SN-117 '+\
            'SN-117M SN-118 SN-119 SN-119M SN-120 SN-121 SN-121M SN-122 '+\
            'SN-123 SN-123M SN-124 SN-125 SN-125M SN-126 SN-127 SN-127M '+\
            'SN-128 SN-128M SN-129 SN-129M SN-130 SN-130M SN-131 SN-131M '+\
            'SN-132 SN-133 SN-134 SN-135 SN-136 SN-137 '
# antimony
nuclides += 'SB-113 SB-115 SB-117 SB-118 SB-118M SB-119 SB-120 SB-120M '+\
            'SB-121 SB-122 SB-122M SB-123 SB-124 SB-124M SB-125 SB-126 '+\
            'SB-126M SB-127 SB-128 SB-128M SB-129 SB-129M SB-130 SB-130M '+\
            'SB-131 SB-132 SB-132M SB-133 SB-134 SB-134M SB-135 SB-136 '+\
            'SB-137 SB-138 SB-139 '
# tellurium
nuclides += 'TE-115 TE-117 TE-118 TE-119 TE-119M TE-120 TE-121 TE-121M '+\
            'TE-122 TE-123 TE-123M TE-124 TE-125 TE-125M TE-126 TE-127 '+\
            'TE-127M TE-128 TE-129 TE-129M TE-130 TE-131 TE-131M TE-132 '+\
            'TE-133 TE-133M TE-134 TE-135 TE-136 TE-137 TE-138 TE-139 '+\
            'TE-140 TE-141 TE-142 '
# iodine
nuclides += 'I-121 I-122 I-123 I-124 I-125 I-126 I-127 I-128 I-129 '+\
            'I-130 I-130M I-131 I-132 I-132M I-133 I-133M I-134 I-134M '+\
            'I-135 I-136 I-136M I-137 I-138 I-139 I-140 I-141 I-142 '+\
            'I-143 I-144 '
# xenon
nuclides += 'XE-122 XE-123 XE-124 XE-125 XE-125M XE-126 XE-127 XE-127M '+\
            'XE-128 XE-129 XE-129M XE-130 XE-131 XE-131M XE-132 XE-133 '+\
            'XE-133M XE-134 XE-134M XE-135 XE-135M XE-136 XE-137 XE-138 '+\
            'XE-139 XE-140 XE-141 XE-142 XE-143 XE-144 XE-145 XE-146 XE-147 '
# cesium
nuclides += 'CS-127 CS-128 CS-129 CS-130 CS-131 CS-132 CS-133 CS-134 '+\
            'CS-134M CS-135 CS-135M CS-136 CS-136M CS-137 CS-138 '+\
            'CS-138M CS-139 CS-140 CS-141 CS-142 CS-143 CS-144 CS-145 '+\
            'CS-146 CS-147 CS-148 CS-149 CS-150 CS-151 '
# barium
nuclides += 'BA-128 BA-129 BA-130 BA-131 BA-131M BA-132 BA-133 BA-133M '+\
            'BA-134 BA-135 BA-135M BA-136 BA-136M BA-137 BA-137M BA-138 '+\
            'BA-139 BA-140 BA-141 BA-142 BA-143 BA-144 BA-145 BA-146 '+\
            'BA-147 BA-148 BA-149 BA-150 BA-151 BA-152 BA-153 '
# lanthanum
nuclides += 'LA-133 LA-134 LA-135 LA-136 LA-137 LA-138 LA-139 LA-140 '+\
            'LA-141 LA-142 LA-143 LA-144 LA-145 LA-146 LA-146M LA-147 '+\
            'LA-148 LA-149 LA-150 LA-151 LA-152 LA-153 LA-154 LA-155 '
# cerium
nuclides += 'CE-134 CE-135 CE-136 CE-137 CE-137M CE-138 CE-139 CE-139M '+\
            'CE-140 CE-141 CE-142 CE-143 CE-144 CE-145 CE-146 CE-147 '+\
            'CE-148 CE-149 CE-150 CE-151 CE-152 CE-153 CE-154 CE-155 '+\
            'CE-156 CE-157 '
# praseodymium
nuclides += 'PR-139 PR-140 PR-141 PR-142 PR-142M PR-143 PR-144 PR-144M '+\
            'PR-145 PR-146 PR-147 PR-148 PR-148M PR-149 PR-150 PR-151 '+\
            'PR-152 PR-153 PR-154 PR-155 PR-156 PR-157 PR-158 PR-159 '
# neodymium
nuclides += 'ND-140 ND-141 ND-141M ND-142 ND-143 ND-144 ND-145 ND-146 '+\
            'ND-147 ND-148 ND-149 ND-150 ND-151 ND-152 ND-153 ND-154 '+\
            'ND-155 ND-156 ND-157 ND-158 ND-159 ND-160 ND-161 '
# promethium
nuclides += 'PM-141 PM-143 PM-144 PM-145 PM-146 PM-147 PM-148 PM-148M '+\
            'PM-149 PM-150 PM-151 PM-152 PM-152M PM-153 PM-154 PM-154M '+\
            'PM-155 PM-156 PM-157 PM-158 PM-159 PM-160 PM-161 PM-162 PM-163 '
# samarium
nuclides += 'SM-143 SM-143M SM-144 SM-145 SM-146 SM-147 SM-148 SM-149 '+\
            'SM-150 SM-151 SM-152 SM-153 SM-154 SM-155 SM-156 SM-157 '+\
            'SM-158 SM-159 SM-160 SM-161 SM-162 SM-163 SM-164 SM-165 '
# europium
nuclides += 'EU-145 EU-146 EU-147 EU-148 EU-149 EU-150 EU-150M EU-151 '+\
            'EU-152 EU-152M EU-153 EU-154 EU-154M EU-155 EU-156 EU-157 '+\
            'EU-158 EU-159 EU-160 EU-161 EU-162 EU-163 EU-164 EU-165 '+\
            'EU-166 EU-167 '
# gadolinium
nuclides += 'GD-146 GD-147 GD-148 GD-149 GD-150 GD-151 GD-152 GD-153 '+\
            'GD-154 GD-155M GD-155 GD-156 GD-157 GD-158 GD-159 GD-160 '+\
            'GD-161 GD-162 GD-163 GD-164 GD-165 GD-166 GD-167 GD-168 GD-169 '
# terbium
nuclides += 'TB-151 TB-152 TB-153 TB-154 TB-154M TB-155 TB-156 TB-156M '+\
            'TB-157 TB-158 TB-158M TB-159 TB-160 TB-161 TB-162 TB-163 '+\
            'TB-164 TB-165 TB-166 TB-167 TB-168 TB-169 TB-170 TB-171 '
# dysprosium
nuclides += 'DY-154 DY-155 DY-156 DY-157 DY-158 DY-159 DY-160 DY-161 '+\
            'DY-162 DY-163 DY-164 DY-165 DY-165M DY-166 DY-167 DY-168 '+\
            'DY-169 DY-170 DY-171 DY-172 '
# holmium
nuclides += 'HO-159 HO-159M HO-160 HO-160M HO-161 HO-161M HO-162 '+\
            'HO-162M HO-163 HO-163M HO-164 HO-164M HO-165 HO-166 '+\
            'HO-166M HO-167 HO-168 HO-169 HO-170 HO-170M HO-171 HO-172 '
# erbium
nuclides += 'ER-160 ER-161 ER-162 ER-163 ER-164 ER-165 ER-166 ER-167 '+\
            'ER-167M ER-168 ER-169 ER-170 ER-171 ER-172 '
# thulium
nuclides += 'TM-165 TM-166 TM-167 TM-168 TM-169 TM-170 TM-171 TM-172 TM-173 '
# ytterbium
nuclides += 'YB-166 YB-167 YB-168 YB-169 YB-169M YB-170 YB-171 YB-172 '+\
            'YB-173 YB-174 YB-175 YB-175M YB-176 YB-177 '
# lutetium
nuclides += 'LU-169 LU-169M LU-170 LU-171 LU-171M LU-172 LU-172M LU-173 '+\
            'LU-174 LU-174M LU-175 LU-176 LU-176M LU-177 LU-177M '
# hafnium
nuclides += 'HF-170 HF-171 HF-172 HF-173 HF-174 HF-175 HF-176 HF-177 '+\
            'HF-177M HF-178 HF-178M HF-179 HF-179M HF-180 HF-180M '+\
            'HF-181 HF-182 '
# tantalum
nuclides += 'TA-177 TA-178 TA-179 TA-180M TA-180 '+\
            'TA-181 TA-182 TA-182M TA-183 '
# tungsten
nuclides += 'W-178 W-180 W-181 W-182 W-183M W-183 W-184 '+\
            'W-185 W-185M W-186 W-187 W-188 W-189 '
# rhenium
nuclides += 'RE-181 RE-182 RE-182M RE-183 RE-184 RE-184M RE-185 '+\
            'RE-186 RE-186M RE-187 RE-188 RE-188M RE-189 '
# osmium
nuclides += 'OS-182 OS-183 OS-184 OS-185 OS-186 OS-187 OS-188 OS-189 '+\
            'OS-189M OS-190 OS-190M OS-191 OS-191M OS-192 OS-193 OS-194 '
# iridium
nuclides += 'IR-185 IR-186 IR-188 IR-189 IR-189M IR-190 IR-191 IR-191M '+\
            'IR-192 IR-192M IR-193 IR-193M IR-194 IR-194M IR-196 IR-196M '
# platinum
nuclides += 'PT-188 PT-189 PT-190 PT-191 PT-192 PT-193 PT-193M PT-194 '+\
            'PT-195 PT-195M PT-196 PT-197 PT-197M PT-198 PT-199 '+\
            'PT-199M PT-200 '
# gold
nuclides += 'AU-193 AU-194 AU-195 AU-195M AU-196 AU-197 AU-197M AU-198 '+\
            'AU-198M AU-199 AU-200 AU-200M '
# mercury
nuclides += 'HG-193 HG-193M HG-194 HG-195 HG-195M HG-196 HG-197 '+\
            'HG-197M HG-198 HG-199 HG-199M HG-200 HG-201 HG-202 '+\
            'HG-203 HG-204 HG-205 HG-206 '
# thallium
nuclides += 'TL-200 TL-201 TL-202 TL-203 TL-204 TL-205 '+\
            'TL-206 TL-207 TL-208 TL-209 TL-210 '
# lead
nuclides += 'PB-200 PB-202 PB-203 PB-204 PB-205 PB-205M PB-206 PB-207 '+\
            'PB-207M PB-208 PB-209 PB-210 PB-211 PB-212 PB-214 '
# bismuth
nuclides += 'BI-205 BI-206 BI-207 BI-208 BI-209 BI-210 BI-210M '+\
            'BI-211 BI-212 BI-212M BI-213 BI-214 '
# polonium
nuclides += 'PO-206 PO-207 PO-208 PO-209 PO-210 PO-211 PO-211M '+\
            'PO-212 PO-213 PO-214 PO-215 PO-216 PO-218 '
# astatine
nuclides += 'AT-216 AT-217 AT-218 '
# radon
nuclides += 'RN-216 RN-217 RN-218 RN-219 RN-220 RN-222 '
# francium
nuclides += 'FR-220 FR-221 FR-222 FR-223 '
# radium
nuclides += 'RA-220 RA-222 RA-223 RA-224 RA-225 RA-226 RA-227 RA-228 '
# actinium
nuclides += 'AC-224 AC-225 AC-226 AC-227 AC-228 '
# thorium
nuclides += 'TH-226 TH-227 TH-228 TH-229 TH-230 TH-231 TH-232 TH-233 TH-234 '
# protactinium
nuclides += 'PA-228 PA-229 PA-230 PA-231 PA-232 '+\
            'PA-233 PA-234M PA-234 PA-235 '
# uranium
nuclides += 'U-230 U-231 U-232 U-233 U-234 U-235 U-235M '+\
            'U-236 U-237 U-238 U-239 U-240 U-241 '
# neptunium
nuclides += 'NP-234 NP-235 NP-236M NP-236 NP-237 '+\
            'NP-238 NP-239 NP-240M NP-240 NP-241 '
# plutonium
nuclides += 'PU-236 PU-237M PU-237 PU-238 PU-239 PU-240 PU-241 '+\
            'PU-242 PU-243 PU-244 PU-245 PU-246 PU-247 '
# americium
nuclides += 'AM-239 AM-240 AM-241 AM-242M AM-242 AM-243 '+\
            'AM-244 AM-244M AM-245 AM-246 AM-246M AM-247 '
# curium
nuclides += 'CM-240 CM-241 CM-242 CM-243 CM-244 CM-245 '+\
            'CM-246 CM-247 CM-248 CM-249 CM-250 CM-251 '
# berkelium
nuclides += 'BK-245 BK-246 BK-247 BK-248 BK-248M BK-249 BK-250 BK-251 '
# californium
nuclides += 'CF-246 CF-248 CF-249 CF-250 CF-251 CF-252 CF-253 CF-254 CF-255 '
# einsteinium
nuclides += 'ES-251 ES-252 ES-253 ES-254M ES-254 ES-255 '
###############################################################################

