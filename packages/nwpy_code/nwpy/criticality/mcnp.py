###############################################################################
"""
    Last edited on June 28, 2019
    
    @author: matz
    
    comments: Prepare the MCNP material input (weight fraction) based
        on stream composition dictionary for mixture of fissile material,
        rock, and water reflected by rock and water
    
    """
###############################################################################
import os
import datetime
import subprocess
from decimal import Decimal
from .. import nuclidedata
###############################################################################


class MCNP(object):
    
    
    """
    A class to write MCNP inputs for nuclear security and safeguards
    material attractiveness FOM evaluation. The MCNP calculations are 
    for bare spheres of nuclear material from different fuel cycles.
    
    """


    @classmethod
    def make_input(cls, filepath, core_mat, core_dens, refl_mat,
                   refl_dens, r_core, kinf=False):
        """Make the MCNP input file
        
        Parameters
        ----------
        filepath: str
            Path to the working directory
        
        core_mat: dict
            Dictionary of the composition of the heavy metal, rock and water
            mixture in the core; isotopes are keys, masses are values
        
        core_dens: float
            Density of the material described in core_mat
        
        refl_mat: dict
            Dictionary of the composition of the rock and water mixture
            in the reflector surrounding the core
        
        refl_dens: float
            Density of the material described in refl_mat
        
        r_core: float
            Radius of the deposition core
            
        kinf (optional): Boolean
            If True, perform a quick kcode calculation with a reflective
            boundary condition around the sphere to determine k-infinity.
            The calculation is made faster by using a kcode entry with 
            fewer neutrons per cycle, fewer total cycles, and fewer 
            discarded cycles.
            
        """
        
        input = cls.make_header()
        input += cls.make_cells(core_dens, refl_dens, kinf)+'\n'
        input += cls.make_surfaces(r_core, kinf)+'\n'
        input += cls.make_data(core_mat, refl_mat, kinf)+'\n'
        input_file = open(filepath, 'w')
        input_file.write(input)
        input_file.close()
        return(filepath)


    @staticmethod
    def make_header():
        """MCNP input file header"""
        
        header =  'C    WRITTEN BY MCNP INPUT GENERATOR'+'\n'
        header += 'C    FAR FIELD DEPOSITION MINIMUM CRITICAL MASS'+'\n'
        header += 'C    USING THE NWPY PYTHON PACKAGE'+'\n'
        header += 'C    '+str(datetime.datetime.today())[:-7]+'\n'
        header += 'C'+'\n'
        header += 'C    MILOS ATZ'+'\n'
        header += 'C    WASTE MANAGEMENT GROUP'+'\n'
        header += 'C    DEPARTMENT OF NUCLEAR ENGINEERING'+'\n'
        header += 'C    UNIVERSITY OF CALIFORNIA, BERKELEY'+'\n'
        header += 'C'+'\n'
        return(header)


    #--------------------------------------------------------------------------
    # CELL CARDS
    #--------------------------------------------------------------------------
    @classmethod
    def make_cells(cls, core_density, refl_density, kinf):
        """Make a single-material, spherical cell cards for MCNP input
        
        Parameters
        ----------
        core_density: float
            Core material density
        
        refl_density: float
            Deposition reflector density
        
        kinf: Boolean
            Indicates whether to calculate for k-infinity or k-effective
            
        """
        
        core_density = '%.4f' % core_density
        refl_density = '%.4f' % refl_density
        cells = 'C    CELL CARDS'+'\n'
        cells +='1    1  -'+core_density+'  -1    imp:n=1'+'\n'
        if(not kinf):
            cells +='2    2  -'+refl_density+'  1 -2  imp:n=1'+'\n'
            cells +='3    0            2     imp:n=0'+'\n'
        return(cells)
        

    #--------------------------------------------------------------------------
    # SURFACE CARDS
    #--------------------------------------------------------------------------
    @classmethod
    def make_surfaces(cls, radius, kinf):
        """Make spherical surface"""
        
        surf = 'C    SURFACE CARDS'+'\n'
        if(kinf): # use big radius to reduce reflected sphere edge effects
            surf += '*1    so  '+'%.6f' % (1000*float(radius)) +'\n'
        else:
            surf += '1    so  '+'%.6f' % radius+'\n'
            surf += '2    so  '+'%.6f' % (radius+100.0)+'\n'
        return(surf)

    
    #--------------------------------------------------------------------------
    # DATA CARDS
    # The MCNP input data cards contain the material information, the
    # criticality control cards
    #--------------------------------------------------------------------------
    @classmethod
    def make_data(cls, m_core, m_refl,  kinf=False):
        """Make the criticality control and material cards
        
        Parameters
        ----------
        particle: str
            Either "n" (neutron) or "p" (photon)
        
        m_core: dict
            Dictionary of the composition of the heavy metal, rock and water
            mixture in the core; isotopes are keys, masses are values
        
        m_refl: dict
            Dictionary of the composition of the rock and water mixture
            in the reflector surrounding the core
            
        kinf (optional): Boolean
            If True, use a kcode entry with fewer neutrons per cycle, fewer
            total cycles, and fewer discarded cycles.
        
        """
        
        data = 'C    DATA CARDS'+'\n'
        data += 'MODE N'+'\n'
        data += cls.make_materials(m_core, 1)+'\n'
        data += cls.make_materials(m_refl, 2)+'\n'
        data += 'C'+'\n'
        # source cards
        if(kinf):
            data += cls.make_kcode(n_per_cycle=250, skip=20, cycles=50)+'\n'
        else:
            data += cls.make_kcode()+'\n'
        return(data)


    #--------------------------------------------------------------------------
    # MATERIAL CARDS
    #--------------------------------------------------------------------------
    @classmethod
    def make_materials(cls, mat_dict, n=1, cutoff=1e-8):
        """Given the composition dictionary from a stream instance, return a
            materials input for MCNP"""
        
        xslib = cls._get_mcnp_data()
        if(n==1):
            mat = 'C    MATERIAL CARDS'+'\n'
        else:
            mat = ''
        mat += 'm'+str(n)+'   '
        mass = sum(mat_dict.values())
        for key in mat_dict.keys():
            zaid = cls._make_zaid(key)
            if(mat_dict[key] > cutoff):
                lib = cls._get_lib(zaid, xslib)
                if(lib is not None):
                    mf = mat_dict[key]/mass
                    mat += zaid+lib+' -'+'%.5E'%Decimal(mat_dict[key]/mass)
                    mat += '\n'+'     '
        mat = mat[:-6] # delete the last five spaces and last new line
        mat += '\n'+'mt'+str(n)+' lwtr.60t' # thermal scattering
        return(mat)


    @staticmethod
    def _make_zaid(nuclide):
        """Given a string that characterizes a radionuclide, make a ZAID"""

        #el, A, meta = nuclidedata.determine_nuclide_info(nuclide)
        el, A, meta = nuclidedata.determine_nuclide_info(nuclide)
        #zaid = str(nuclidedata.Z[el])
        zaid = str(Z[el])
        if(len(A) < 2):
            zaid += '0'
        if(len(A) < 3):
            zaid += '0'
        zaid += str(A)
        return(zaid)
    
    
    @staticmethod
    def _get_lib(zaid, obj):
        """Check if xs exists with good lib
        
        Parameters
        ----------
        zaid: str
            Nuclide Z-A-ID
        
        obj: str OR dict
            If str: text from MCNP_DATA/xslib
            If dict: from xslib module
            
        """
        
        lib = ['.80c', '.70c', '.66c', '.55c', '.42c', '.24c']
        found = False
        i = 0
        while(found==False and i < len(lib)-1):
            id = zaid+lib[i]
            if(id in obj):
                found = True
                return(lib[i])
            else:
                i += 1
        

    @staticmethod
    def _get_mcnp_data():
        """Get the path to the MCNP data"""
        
        path = os.environ['DATAPATH'].split(os.pathsep)
        path = [x for x in path if 'MCNP' in x]
        #path = [x for x in path if 'DATA' in x]
        return(open(os.path.join(path[0], 'xsdir')).read()) # expect only one entry


    #--------------------------------------------------------------------------
    # SOURCE CARDS
    #--------------------------------------------------------------------------
    @staticmethod
    def make_kcode(n_per_cycle=5000, skip=50, cycles=100):
        """Make some friggin neutrons (kcode, ksrc cards)"""

        neutrons = 'C    CRITICALITY CONTROL CARD'+'\n'
        neutrons += ('kcode'+' '+str(n_per_cycle)+
                     ' 1.0 '+str(skip)+' '+str(cycles)+'\n')
        neutrons += 'ksrc  0 0 0' # put 'em in the middle
        return(neutrons)


    ###########################################################################
    # RUNNING MCNP
    # Parallelization?
    ###########################################################################
    @classmethod
    def run_mcnp(cls, path_to_infile, print_to_terminal=False, version=6):
        """Run MCNP on the requested input file"""
        
        path_to_mcnp = cls.find_mcnp(version)
        if(print_to_terminal):
            subprocess.call([path_to_mcnp, 'name='+path_to_infile])
        else:
            devnull = open(os.devnull, 'w')
            subprocess.call([path_to_mcnp, 'name='+path_to_infile],
                            stdout=devnull, stderr=devnull)
        return(path_to_infile+'o') # outfile


    @classmethod
    def parallel_run_mcnp(cls, path_to_infile, nodes=1, cores=20):
        """Run MCNP6 in parallel on the requested input file; requires this
        Python script (serial) to have been initialized on parallel resources
        
        This is intended to run on Savio, where there are 20 cores per node.
        
        """
        
        path_to_mcnp = cls.find_mcnp(version=6)+'.mpi'
        devnull = open(os.devnull, 'w')
        subprocess.call(['mpirun', path_to_mcnp, 'name='+path_to_infile,
                         'tasks', str(nodes*cores)], stdout=devnull,
                        stderr=devnull)
        return(path_to_infile+'o') # outfile


    @staticmethod
    def find_mcnp(version=6):
        """Find the path to the mcnp executable based on the contents
        of the .bashrc path environment variables"""
        
        mcnp_exec = 'mcnp'+str(version)
        path_to_mcnp = []
        for p in os.environ['PATH'].split(os.pathsep):
            try: # if not a directory, this will fail
                temp = os.listdir(p)
            except OSError:
                continue
            else:
                if(mcnp_exec in os.listdir(p)): # unix
                    path_to_mcnp.append(os.path.join(p, mcnp_exec))
                elif(mcnp_exec+'.exe' in os.listdir(p)): # windows
                    path_to_scale.append(os.path.join(p, mcnp_exec+'.exe'))
                else:
                    continue
        warning = ('Path to scalerte executable not found. Make sure '+
                   'it is added to .bashrc PATH environment')
        assert path_to_mcnp != [], warning
        assert os.path.exists(path_to_mcnp[0]), warning
        return(path_to_mcnp[0])


###############################################################################
# NUCLIDE IDENTIFIERS (from nwpy.fuelcycle.nuclidedata)
# Parse nuclide ID of the form EEAAM, where "EE" is the element (single "E"
# okay); "AA" is mass (single "A" okay); "M" indicates isotope metastability.
###############################################################################


Z = {'h':1, 'he':2, 'li':3, 'be':4, 'b':5, 'c':6, 'n':7, 'o':8, 'f':9,
     'ne':10, 'na':11, 'mg':12, 'al':13, 'si':14, 'p':15, 's':16, 'cl':17,
     'ar':18, 'k':19, 'ca':20, 'sc':21, 'ti':22, 'v':23, 'cr':24, 'mn':25,
     'fe':26, 'co':27, 'ni':28, 'cu':29, 'zn':30, 'ga':31, 'ge':32, 'as':33,
	 'se':34, 'br':35, 'kr':36, 'rb':37, 'sr':38, 'y':39, 'zr':40, 'nb':41,
     'mo':42, 'tc':43, 'ru':44, 'rh':45, 'pd':46, 'ag':47, 'cd':48, 'in':49,
     'sn':50, 'sb':51, 'te':52, 'i':53, 'xe':54, 'cs':55, 'ba':56, 'la':57,
     'ce':58, 'pr':59, 'nd':60, 'pm':61, 'sm':62, 'eu':63, 'gd':64, 'tb':65,
     'dy':66, 'ho':67, 'er':68, 'tm':69, 'yb':70, 'lu':71, 'hf':72, 'ta':73,
	 'w':74, 're':75, 'os':76, 'ir':77, 'pt':78, 'au':79, 'hg':80, 'tl':81,
     'pb':82, 'bi':83, 'po':84, 'at':85, 'rn':86, 'fr':87, 'ra':88, 'ac':89,
     'th':90, 'pa':91, 'u':92, 'np':93, 'pu':94, 'am':95, 'cm':96, 'bk':97,
     'cf':98, 'es':99, 'fm':100, 'md':101, 'no':102, 'lr':103}


###############################################################################
