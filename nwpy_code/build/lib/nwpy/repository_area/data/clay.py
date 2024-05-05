###################################################################################
"""
    Last edited on September 15, 2018
    
    @author: matz
    
    comments: Data for the calculation of heat in a clay repository
              as per the SNL/LLNL 2011 report.
    
    """
###################################################################################
# CLAY REPOSITORY
###################################################################################
# Thermal constraints
constraint = {}
constraint['wp'] = 100.0     # degrees C at the surface of the WP
constraint['rock'] = 100.0   # degrees C in rock
###################################################################################
# Repository Design
# In SNL/LLNL 2011, both the SNF and HLW design concepts have drifts spaced
# 20 m apart. However, in the SNF design, the packages are spaced 10 m apart,
# while in the HLW design, the packages are spaced 6 m apart. Because this code
# alters the spacing between drifts and packages, a single design was chosen to
# be that of the SNF.
spacing = {}
spacing['drift'] = 30.0 # m
spacing['pkg'] = 10.0 # m
# EBS dimensions for spent nuclear fuel
ebs = {}
ebs['snf'] = {}
ebs['snf']['layer'] = ['liner', 'envelope', 'buffer', 'overpack']
ebs['snf']['material'] = ['steel', 'carbon steel', 'bentonite', 'carbon steel']
ebs['snf']['dr'] = [0.025, 0.006, 0.8, 0.05]
# SNF package is 0.82 m diameter
# SNF canister is 0.72 m diameter
# SNF overpack thickness is 0.05 m
# EBS dimensions for high level waste
ebs['hlw'] = {}
ebs['hlw']['layer'] = ['liner', 'overpack']
ebs['hlw']['material'] = ['steel', 'carbon steel']
ebs['hlw']['dr'] = [0.01, 0.055]
# HLW canister is 0.61 m diameter
# HLW package is 0.72 m diameter
# Overpack thickness = (0.72-0.61)/2 = 0.055 m
###################################################################################