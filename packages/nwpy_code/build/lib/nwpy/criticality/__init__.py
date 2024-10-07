###############################################################################
"""
    Last edited on June 20, 2019
    
    @author: matz
    
    comments: Package initialization
    
"""
###############################################################################


def renormalize(comp):
    """Renormalize a dictionary so that the values sum to one
    
    Parameters
    ----------
    comp: dict
        Dictionary with values to be renormalized
    
    Returns
    -------
    Dictionary with values normalized to 1.0
    
    """
        
    total = sum(comp.values())
    return(dict((k,v/total) for k,v in comp.iteritems()))



###############################################################################
