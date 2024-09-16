###################################################################################
"""
    Last edited on September 15, 2018
	
    @author: matz
	
    comments: Data for the calculation of heat in a granite repository 
              as per the SNL/LLNL 2011 report.
	
"""
###################################################################################
# GRANITE REPOSITORY
###################################################################################
# Thermal constraints
constraint = {}
constraint['wp'] = 100.0     # degrees C at the surface of the WP
constraint['rock'] = 200.0   # degrees C in rock
###################################################################################
# Repository Design
# In SNL/LLNL 2011, both the SNF and HLW design concepts have drifts spaced 20 m
# apart and the packages spaced 10 m apart.
spacing = {}
spacing['drift'] = 20.0 # m
spacing['pkg'] = 10.0 # m
# EBS dimensions for spent nuclear fuel
ebs = {}
ebs['snf'] = {}
ebs['snf']['layer'] = ['buffer', 'overpack']
ebs['snf']['material'] = ['bentonite', 'copper']
ebs['snf']['dr'] = [0.35, 0.05]
# SNF package is 0.82 m diameter
# SNF canister is 0.72 m diameter
# SNF overpack thickness is 0.05 m
# EBS dimensions for high level waste
ebs['hlw'] = {}
ebs['hlw']['layer'] = ['buffer', 'overpack']
ebs['hlw']['material'] = ['bentonite', 'copper']
ebs['hlw']['dr'] = [0.35, 0.105]
# HLW canister diameter is 0.61 m
# HLW package diameter is 0.82 m
# HLW overpack thickness is 0.105 m
###################################################################################
