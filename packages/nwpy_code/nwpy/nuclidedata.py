###############################################################################
"""
    Last modified on May 13, 2019

    @author: matz

    comment: data file containing actinide list and element z-values

"""
###############################################################################
# Define element and nuclide species groups
actinides = ['ac', 'th', 'pa', 'u', 'np', 'pu', 'am', 'cm',
             'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr']
minor_actinides = ['np', 'am', 'cm', 'bk', 'cf', 'es',
                   'fm', 'md', 'no', 'lr']
transuranics = ['np', 'pu', 'am', 'cm', 'bk', 'cf',
                'es', 'fm', 'md', 'no', 'lr']
u3 = ['pa', 'u'] # denatured/impure u233 bred from th
llfp = ['se79', 'tc99', 'sn126', 'i129', 'cs135']
###############################################################################
# We need to be able to build a library of ZAIDs for the imported list of
# nuclides. We'll have the A value and the atomic abbreviation, but we need
# to know the Z for each element.
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
# Parse nuclide ID of the form EEAAM, where "EE" is the element (single "E"
# okay); "AA" is mass (single "A" okay); "M" indicates isotope metastability.
def determine_nuclide_info(nuclide_id):
    """Parse information contained in nuclide ID
        
    Parameters
    ----------
    nuclide_id: str
        Isotopic nuclide identifier of the form E(E)AA(A)M or EE-AAAM
    
    Returns
    -------
    Tuple with three str containing
        1. Element symbol
        2. Mass value
        3. Indication whether the nuclide is a metastable state
    
    """
    
    assert len(nuclide_id)<=6, ('Input must be str and should be '+\
                                'no longer than 6 characters')
    nuclide_id = nuclide_id.lower() # using lowercase letters
    el = []
    mass = []
    m = 0
    meta=False
    for i in range(0, len(nuclide_id)):
        char = nuclide_id[i]
        if(char.isalpha() and i < 2):
            el.append(char)
        elif(char.isdigit()):
            mass.append(char)
        elif(char.lower()=='m' and i>=2):
            meta=True
        else:
            continue
    if(meta==True):
        return(''.join(el), ''.join(mass), 'meta')
    else:
        return(''.join(el), ''.join(mass), 'not meta')


###############################################################################
# Functions to assess nuclide grouping
###############################################################################
def is_group(nuclide, group):
    """Determine if nuclide is a member of a specific group
        
    Parameters
    ----------
    nuclide: str
        Nuclide identifier of the form E(E)AA(A)M
    
    group: str
        Group or species for which to check if nuclide is a member
        
    Returns
    -------
    Boolean indicating whether or not the nuclide is a member of the group
    
    """
    
    group = group.lower()
    if(any(x in group for x in ['fp', 'fissionproduct', 'fission product'])):
        if(any(x in group for x in ['llfp', 'long-lived', 'long lived'])):
            return(is_llfp(nuclide))
        else:
            return(is_fissionproduct(nuclide))
    elif(group=='act' or group=='actinide'):
        return(is_actinide(nuclide))
    elif(group=='tru' or 'transuranic' in group):
        return(is_transuranic(nuclide))
    elif(group=='ma' or 'minor' and 'actinide' in group):
        return(is_minoractinide(nuclide))
    else:
        el, A, meta = determine_nuclide_info(nuclide)
        if(group=='u3'):
            if(el in ['u', 'pa']):
                return(True)
        else:
            if(el==group):
                return(True)
            else:
                return(False)


def is_fissionproduct(nuclide):
    """Determine if nuclide is a fission product
    
    Parameters
    ----------
    nuclide: str
        Nuclide identifier of the form E(E)AA(A)M
        
    Returns
    -------
    Boolean
    
    """
    
    # lightest/heaviest = V/HF (Z=23, 72); ENDF
    el, A, meta = determine_nuclide_info(nuclide)
    if(Z[el] <= Z['hf'] and Z[el] >= Z['v']):
        return(True)
    else:
        return(False)


def is_actinide(nuclide):
    """Determine if nuclide is an actinide
        
    Parameters
    ----------
    nuclide: str
        Nuclide identifier of the form E(E)AA(A)M
    
    Returns
    -------
    Boolean
    
    """

    el, A, meta = determine_nuclide_info(nuclide)
    if(el in actinides):
        return(True)
    else:
        return(False)


def is_transuranic(nuclide):
    """Determine if nuclide is a transuranic

    Parameters
    ----------
    nuclide: str
        Nuclide identifier of the form E(E)AA(A)M
        
    Returns
    -------
    Boolean
    
    """
    
    el, A, meta = determine_nuclide_info(nuclide)
    if(el in transuranics):
        return(True)
    else:
        return(False)


def is_minoractinide(nuclide):
    """Determine if nuclide is a minor actinide
        
    Parameters
    ----------
    nuclide: str
        Nuclide identifier of the form E(E)AA(A)M
    
    Returns
    -------
    Boolean
        
    """
    
    el, A, meta = determine_nuclide_info(nuclide)
    if(el in minor_actinides):
        return(True)
    else:
        return(False)


def is_llfp(nuclide):
    """Determine if nuclide is a long-lived fission product
    
    Parameters
    ----------
    nuclide: str
        Nuclide identifier of the form E(E)AA(A)M
        
    Returns
    -------
    Boolean
    
    """
    if(nuclide in llfp):
        return(True)
    else:
        return(False)


###############################################################################


def group_daughters(srs):
    """Group nuclides unsupported by ORIGEN with their decay daughters;
    this method is used in the Reactor class when importing isotopic data
    from CSVs and in the Origen class when writing input files."""
        
    to_drop = []
    for j in range(0, len(srs)):
        nuc = srs.index[j].lower()
        mass = srs[srs.index[j]]
        if(nuc in group_nuclides.keys()):
            to_drop.append(nuc)
            for daughter in group_nuclides[nuc].keys():
                branchfrac = group_nuclides[nuc][daughter]
                try: # assume daughter is already in composition df
                    srs[daughter] += mass*branchfrac
                except KeyError: # if daughter not in df, make new entry
                    srs[daughter] = mass*branchfrac
    srs = srs.drop(to_drop)
    return(srs)


# No data in ORIGEN decay library for these nuclides - lump w/ decay daughters
group_nuclides = {}
# Cu-81 decays by (1) beta decay to Zn-81, and (2) beta+n decay to Zn-80
group_nuclides['cu81'] = {'zn81': 0.4705, 'zn80': 0.5295}
# Se-85m decays by beta decay to Br-85
group_nuclides['se85m'] = {'br85': 1.0}
# Br-86m decays by beta decay to Kr-86
group_nuclides['br86m'] = {'kr86': 1.0}
# Rh-109m decays by isomeric transition to Rh-109
group_nuclides['rh109m'] = {'rh109': 1.0}
# Rh-123 decays by (1) beta decay to Pd-123 and (2) beta+n decay to Pd-122
group_nuclides['rh123'] = {'pd123': 0.8289, 'pd122': 0.1711}
# Pd-125 decays by beta decay to Ag-125
group_nuclides['pd125'] = {'ag125': 1.0}
# Pd-126 decays by (1) beta decay to Ag-126, and (2) beta+n decay to Ag-125
group_nuclides['pd126'] = {'ag126': 0.9497, 'ag125': 0.0503}
# I-145 decays by beta decay to Xe-145
group_nuclides['i145'] = {'xe145': 1.0}
# Gd-153m decays to Gd-153 by gamma decay
group_nuclides['gd153m'] = {'gd153': 1.0}
# Tb-162m decays by isomeric transition to Tb-162
# http://www.radiochemistry.org/periodictable/elements/isotopes_data/65.html
group_nuclides['tb162m'] = {'tb162': 1.0}
# Tb-163m decays by isomeric transition to Tb-163
group_nuclides['tb163m'] = {'tb163': 1.0}
# Tm-170m decay by isomeric transition to Tm-170
group_nuclides['tm170m'] = {'tm170': 1.0}
# Pt-186 decays by beta+ decay to Ir-186
group_nuclides['pt186'] = {'ir186': 1.0}
# Pt-187 decays by beta+ decay to Ir-187
group_nuclides['pt187'] = {'ir187': 1.0}
# Tl-196 decays by beta+ decay to Hg-197
group_nuclides['tl196'] = {'hg196': 1.0}
# Tl-197 decays by beta+ decay to Hg-197
group_nuclides['tl197'] = {'hg197': 1.0}
# Pa-227 decays by (1) alpha to Ac-223, and (2) isomeric transition to Th-227
group_nuclides['pa227'] = {'ac223':0.85, 'th227':0.15}
# SCALE ON SAVIO DOES NOT LIKE F-19 AND IT SHOULDN'T MATTER
group_nuclides['f19'] = {'f19':0.0} # make it disappear
# Known others:
# tl194
# nd139m
# sb116m
# tb150
# tl194m
# ta172
# ta173
# ta174
# ta175
# ir187
# ir184


###############################################################################
