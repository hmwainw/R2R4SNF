###############################################################################
"""
    Last edited on June 28, 2019
    
    @author: matz
    
    comments: Process water material inputs for critical mass evaluation
    of mixture of rock, water, and fissile heavy metal.
    
"""
###############################################################################
from rock import Rock
###############################################################################


class Water(Rock):
    
    
    """
    The water in the geologic formation in which the fissile material 
    precipitate exists and for which the critical mass is being evaluated. 
    This object calculates the masses of each element and oxygen in the 
    pore water.
    
    Although the methods here are similar to those in the HeavyMetal object,
    they are different in that the composition of the rock is specified by 
    the mineral oxides, which need to be parsed somewhat differently. The
    data here is in the same form as for the Rock object so they share
    methods.
    
    """
    
    
    def __init__(self):
        """Initialize heavy metal material instance"""
        
        self.comp = water
        self.density = density
        self.elemental_mf()
        self.elemental_mm()
    

###############################################################################
# Water composition
water = {}
water['h2o'] = 1.0000
#water['cl'] = 0.0000
###############################################################################
# Water density
density = 1.0 # g/cc
###############################################################################