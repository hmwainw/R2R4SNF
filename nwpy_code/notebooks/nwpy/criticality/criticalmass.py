###############################################################################
"""
    Last edited on June 28, 2019
    
    @author: matz
    
    comments: Calculate the minimum critical mass for depositions of fissile
        material in the far field of a geologic repository for nuclear waste
    
"""
###############################################################################
import numpy as np
from collections import Counter
import rock
import water
import heavymetal
import __init__
import data
from criticalvolume import CriticalVolume
from .. import nuclidedata
###############################################################################


class CriticalMass(object):


    """
    This class performs a parametric sweep over the heavy metal volume
    fraction to determine minimum critical mass for some input heavy
    metal composition.

    """


    @classmethod
    def calculate_minimum(cls, comp, porosity=0.1, mobility={'u':1.0},
                          dhmvf=0.02, name='', wdpath='./'):
        """Calculate the minimum critical mass for a given heavy metal 
        composition and host rock porosity.
        
        Paramters
        ---------
        comp: dict or nwpy Stream instance
            Data for heavy metal material composition (by mass)
            
        porosity (optional): float
            Void volume in the deposition host rock that can be occupied 
            by heavy metal and/or water
        
        mobility (optional): dict
            Dictionary describing the relative mobility of heavy metal
            elements to one another - informs the extent to which elements
            in the input heavy metal stream are included in the deposition.
            The default is that uranium is the only mobile species; if, for
            example, you want to include Pu to the same extent as U, you
            would use: mobility={'u':1.0, 'pu':1.0}. If you want to account
            for the fact that Pu is less mobile than U (say, by an order of
            magnitude), you might say: mobility={'u':1.0, 'pu':0.1}.
        
        name (optional): str
            Base filename for MCNP file I/O
        
        wdpath (optional): str
            Path to working directory for MCNP file I/O
            
        Returns
        -------
        Minimum critical mass of a deposition of the requested heavy metal,
        rock (sandstone), and water.
        
        """
        
        cls.cm = []
        cls.hmvf = []
        comp = cls.mobility_adjustment(comp, mobility)
        hm = heavymetal.HeavyMetal(comp)
        r = rock.Rock()
        w = water.Water()
        hmvf = np.arange(dhmvf, porosity+dhmvf, dhmvf)
        if(len(name)>0):
            name += '_'
        for i in xrange(0, len(hmvf)):
            x = hmvf[i]
            casename = name+'hmvf'+str(i)
            mass = cls.calculate_cm(casename, hm, r, w, porosity, x, wdpath)
            if(mass is not None):
                cls.cm.append(mass)
                cls.hmvf.append(x)
            if(not cls.continue_search(cls.cm)):
                break
        if(len(cls.cm) > 0):
            return(min(cls.cm))
        else: # criticality screened
            return


    @staticmethod
    def continue_search(cm, buffer=3):
        """Review the results from the previous iterations to end the
        sweep over VVF when the minimum is found. Assumes the critical
        mass will decrease with increasing HMVF from HMVF=0 until a 
        minimum is found.
        
        Parameters
        ----------
        cm: list
            Critical masses accumulated in the parametric sweep
            
        buffer (optional): int
            The number of cases past the minimum to include to ensure 
            that the minimum is correctly captured.
        
        Returns
        -------
        Boolean indicating whether to continue the sweep
        
        """
        
        if(len(cm) < buffer+1): # list is too short, need more values
            return(True)
        else: # make sure minimum is sufficiently buffered
            test = []
            for i in xrange(1, buffer+1):
                test.append(min(cm)!=cm[-i])
            if(all(test)):
                return(False)
            else:
                return(True)
    

    @staticmethod
    def mobility_adjustment(hmcomp, mobility):
        """Adjust the composition of heavy metal based on user inputs for
        relative mobility of actinide elements"""
        
        hm_dep = {}
        for nuc in hmcomp:
            el = nuclidedata.determine_nuclide_info(nuc)[0]
            if(el in mobility.keys()):
                hm_dep[nuc] = hmcomp[nuc]*mobility[el]
            else:
                continue
        return(hm_dep) #__init__.renormalize(hm_dep))


    @classmethod
    def calculate_cm(cls, name, hm, r, w, porosity, hmvf, wdpath='./'):
        """Calculate the critical mass of a single deposition configuration"""
        
        core = cls.make_core(hm, r, w, porosity, hmvf)
        ref = cls.make_refl(r, w, porosity)
        cv=CriticalVolume.calculate(name,core[0],core[1],ref[0],ref[1],wdpath)
        if(cv is not None):
            return(cv*hmvf*hm.density) # grams
        else:
            return(None)


    @staticmethod
    def make_core(hm, r, w, vvf, hmvf):
        """Generate the composition and density of the deposition core region
        given the compositions, densities, and volume fractions of the rock, 
        water, and heavy metal"""
        
        r_mat, r_dens = r.make(1.0-vvf)
        w_mat, w_dens = w.make(1.0-vvf-hmvf)
        hm_mat, hm_dens = hm.make(hmvf)
        density = r_dens*(1.0-vvf)+w_dens*(1.0-vvf-hmvf)+hm_dens*hmvf
        return(dict(Counter(r_mat)+Counter(w_mat)+Counter(hm_mat)), density)
    
    
    @staticmethod
    def make_refl(r, w, vvf):
        """Generate the composition and density of the deposition reflector
        region given the compositions, densities, and volume fractions of the 
        rock and water"""
        
        r_mat, r_dens = r.make(1.0-vvf)
        w_mat, w_dens = w.make(vvf)
        density = r_dens*(1.0-vvf)+w_dens*(1.0-vvf)
        return(dict(Counter(r_mat)+Counter(w_mat)), density)


###############################################################################
#import rock
#import water
#import heavymetal
#comp = {'u235':5.0, 'u238':95.0}
#hm = heavymetal.HeavyMetal(comp)
#r=rock.Rock()
#w = water.Water()
#criticalmass.CriticalMass.calculate_cm('test', hm, r, w, 0.8, 0.12, wdpath='./fake/')
