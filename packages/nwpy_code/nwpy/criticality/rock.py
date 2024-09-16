###############################################################################
"""
    Last edited on June 17, 2019
    
    @author: matz
    
    comments: Process rock material inputs for critical mass evaluation 
    of mixture of rock, water, and fissile heavy metal.
    
"""
###############################################################################
import numpy as np
import data
from .. import nuclidedata
v1_radius = (0.75/np.pi)**(1/3.0) # radius (cm) for sphere with V = 1 (cc)
###############################################################################


class Rock(object):
    
    
    """
    The rock in which the fissile material precipitate exists and for which
    the critical mass is being evaluated. This object calculates the masses
    of each element and oxygen in the deposition.
    
    Although the methods here are similar to those in the HeavyMetal object,
    they are different in that the composition of the rock is specified by 
    the mineral oxides, which need to be parsed somewhat differently.
    
    """
    
    
    def __init__(self):
        """Initialize heavy metal material instance"""
    
        self.comp = rock
        self.density = density
        self.elemental_mf()
        self.elemental_mm()
    
    
    def make(self, vf, radius=v1_radius):
        """Get the rock isotopes that will go into the deposition
        based on the deposition radius and rock volume fraction
        
        Parameters
        ----------
        vf: float
            Value between 0 and 1 characterizing the volume 
            fraction of rock in the deposition.
        
        radius (optional): float
            Radius of the deposition
        
        Returns
        -------
        1. Dictionary of rock isotope masses in the precipitate (incl oxygen)
        2. Density of the rock mineral (oxide) precipitate
        
        """
        
        v_dep = (4/3.0)*np.pi*radius**3
        mat = {}
        mat['o16'] = 0.0
        for oxide in self.comp.keys():
            stoich = self._get_stoichiometry(oxide)
            el = self._get_element(oxide)
            mm_ox = stoich[0]*self.mm[el]+stoich[1]*15.999
            for iso in self.mf[el].keys():
                if(not iso in mat.keys()):
                    mat[iso] = 0.0
                mat[iso] += (self.mf[el][iso]*self.mm[el]*stoich[0]*
                             self.comp[oxide]*self.density*vf*v_dep/mm_ox)
            mat['o16'] += (v_dep*vf*self.density*self.comp[oxide]*
                           stoich[1]*15.999/mm_ox)
        return(mat, density)
           
        
    def elemental_mf(self):
        """Get isotopic fractionation of elements"""
    
        self.mf = {}
        for ox in self.comp.keys():
            el = self._get_element(ox)
            if(el not in self.mf):
                self.mf[el] = {}
            for iso in data.mm.keys():
                eli = nuclidedata.determine_nuclide_info(iso)[0]
                if(eli == el):
                    self.mf[el][iso] = data.a[iso]
    
    
    def elemental_mm(self):
        """Calculate the elemental molar masses of each element in the rock"""
    
        self.mm = {}
        for oxide in self.comp.keys():
            el = self._get_element(oxide)
            if(el not in self.mm.keys()):
                self.mm[el] = 0.0
                for iso in data.mm.keys():
                    eli = nuclidedata.determine_nuclide_info(iso)[0]
                    if(eli == el):
                        self.mm[el] += data.a[iso]/data.mm[iso]
        for el in self.mm.keys():
            self.mm[el] = 1/self.mm[el]


    @staticmethod
    def _get_element(oxide):
        """Given a string characterizing a mineral oxide, 
        return the non-oxygen element in the mineral"""
        
        el = ''
        for s in oxide:
            if(s is not 'o'):
                try:
                    float(s)
                    break
                except ValueError:
                    el += s
        return(el)
        
    
    @staticmethod
    def _get_stoichiometry(oxide):
        """Given a string characterizing a mineral oxide, return
        the stoichiometry of the two elements that comprise it"""
        
        stoich = []
        tag = oxide.find('o')
        mineral = oxide[0:tag]
        try:
            stoich.append(int(mineral[-1]))
        except ValueError:
            stoich.append(1)
        oxygen = oxide[tag:]
        try:
            stoich.append(int(oxygen[-1]))
        except ValueError:
            stoich.append(1)
        return(stoich)
        

###############################################################################
# Rock composition
rock = {}
rock['sio2']  = 0.7870
rock['tio2']  = 0.0025
rock['al2o3'] = 0.0480
rock['fe2o3'] = 0.0110
rock['feo']   = 0.0030
rock['mno']   = 0.0003
rock['mgo']   = 0.0120
rock['cao']   = 0.0550
rock['na2o']  = 0.0045
rock['k2o']   = 0.0130
rock['h2o']   = 0.0130
rock['p2o5']  = 0.0008
rock['co2']   = 0.0500
rock['so3']   = 0.0007
###############################################################################
# Rock grain density
density = 2.71 # g/cc
###############################################################################
