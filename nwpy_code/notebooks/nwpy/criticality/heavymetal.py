###############################################################################
"""
    Last edited on June 18, 2019
    
    @author: matz
    
    comments: Process heavy metal material inputs for critical mass 
    evaluation of mixture of rock, water, and fissile heavy metal.
    
"""
###############################################################################
import numpy as np
import data
from .. import nuclidedata
v1_radius = (0.75/np.pi)**(1/3.0) # radius (cm) for sphere with V = 1 (cc)
###############################################################################


class HeavyMetal(object):


    """
    Heavy metal precipitate in the far-field deposition based on isotopic 
    composition input by user. This object calculates the elemental and
    oxide molar masses of the heavy metal elements and based on the 
    desired deposition radius and volume fraction of heavy metal returns
    the masses of each element and oxide in the deposition
    
    """


    def __init__(self, composition):
        """Initialize heavy metal material instance
        
        Parameters
        ----------
        composition: dict
            Isotopic composition by mass

        """
        
        self.name = 'hm'
        self.comp = composition
        # get mixture attributes
        self.elemental_mf() # elemental isotopic fractionation
        self.elemental_mm() # elemental molar mass
        self.oxide_mm() # oxide molar masses
        self.oxide_mf() # oxide mass fractions
        self.get_density() # oxide mixture density
    
    
    def make(self, hmvf, radius=v1_radius):
        """Get the heavy metal materials that will go into the deposition
        based on the deposition radius and heavy metal volume fraction
        
        Parameters
        ----------
        hmvf: float
            Value between 0 and 1 characterizing the volume fraction of
            heavy metal in the deposition.
        
        radius (optional): float
            Radius of the deposition
        
        Returns
        -------
        1. Dictionary of heavy metal masses in the precipitate (incl oxygen)
        2. Density of the heavy metal oxide precipitate
        
        """
        
        v_dep = (4/3.0)*np.pi*radius**3
        mat = {}
        mat['o16'] = 0.0
        for el in self.mf.keys():
            for iso in self.mf[el].keys():
                mat[iso] = (self.mf[el][iso]*self.mm[el]*oxide[el][0]*
                            self.mf_ox[el]*self.density*hmvf*v_dep/
                            self.mm_ox[el])
            mat['o16'] += (v_dep*hmvf*self.density*self.mf_ox[el]*
                           oxide[el][1]*15.999/self.mm_ox[el])
        return(mat, self.density)
            
    
    def get_density(self):
        """Calculate density of the mixture of oxides based 
        on their mass fractions in the precipitate"""
        
        self.density = 0.0
        for el in self.mf_ox.keys():
            self.density += self.mf_ox[el]*density[el]


    def elemental_mf(self):
        """Calculate the isotopic fractionation of each element"""
    
        self.mf = {}
        for iso in self.comp.keys(): # iso is a nuclide string i.e. "u235"
            el = nuclidedata.determine_nuclide_info(iso)[0]
            if(el not in self.mf.keys()):
                self.mf[el] = {}
            self.mf[el][iso] = self.comp[iso]
        for el in self.mf.keys():
            tot = sum(self.mf[el].values())
            for iso in self.mf[el].keys():
                self.mf[el][iso] = self.mf[el][iso]/tot
    
    
    def elemental_mm(self):
        """Calculate elemental molar mass for each heavy metal element"""
        
        self.mm = {}
        moles = {}
        masses = {}
        for iso in self.comp.keys(): # iso is a nuclide string i.e. "u235"
            el = nuclidedata.determine_nuclide_info(iso)[0]
            if(el not in moles.keys()):
                moles[el] = 0.0
                masses[el] = 0.0
            masses[el] += self.comp[iso]
            moles[el] += self.comp[iso]/data.mm[iso]
        for el in moles.keys():
            self.mm[el] = masses[el]/moles[el]


    def oxide_mm(self):
        """Calculate the molar mass of the oxide of each element"""
        
        self.mm_ox = {}
        for el in self.mm:
            self.mm_ox[el] = oxide[el][0]*self.mm[el] + oxide[el][1]*15.999


    def oxide_mf(self):
        """Calculate the mass fractions of each oxide 
        with respect to the total oxide mixture"""
        
        self.mf_ox = {}
        temp = 0.0
        for el in self.mm.keys():
            self.mf_ox[el] = 0.0
            for iso in self.comp:
                eli = nuclidedata.determine_nuclide_info(iso)[0]
                if(eli == el):
                    self.mf_ox[el]+= (self.comp[iso]*self.mm_ox[el]/self.mm[el]/oxide[el][0])
        m_oxide = sum(self.mf_ox.values())
        for el in self.mf_ox.keys():
            self.mf_ox[el] = self.mf_ox[el]/m_oxide


###############################################################################
# Oxide speciations
oxide = {}
oxide['u']  = [1, 2] # UO2
oxide['np'] = [1, 2] # NpO2
oxide['pu'] = [1, 2] # PuO2
oxide['am'] = [2, 3] # Am2O3
oxide['cm'] = [2, 3] # Cm2O3
###############################################################################
# Oxide densities
density = {}
density['u'] =  10.97 # g/cc UO2
density['np'] = 11.10 # g/cc NpO2
density['pu'] = 11.50 # g/cc PuO2
density['am'] = 11.77 # g/cc Am2O3 [1]
density['cm'] = 10.80 # g/cc Cm2O3 [2]
# [1] https://www.atsdr.cdc.gov/toxprofiles/tp156-c4.pdf
# [2] https://webwiser.nlm.nih.gov/getSubstanceData.do?substanceId=413&displaySubstanceName=Curium%20Radioactive&STCCID=&UNNAID=&selectedDataMenuItemID=44&catId=51
###############################################################################





