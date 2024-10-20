###################################################################################
"""
    Last edited on September 15, 2018
    
    @author: matz
    
    comments: Data for the calculation of heat in a salt repository
    as per the SNL/LLNL 2011 report.
    
"""
###################################################################################
# SALT REPOSITORY
###################################################################################
# Thermal constraints
constraint = {}
constraint['wp'] = 200.0     # degrees C at the surface of the WP
constraint['rock'] = 200.0   # degrees C in rock
###################################################################################
# Repository Design
# In SNL/LLNL 2011, both the SNF and HLW design concepts have drifts spaced 20 m
# apart and the packages spaced 20 m apart.
spacing = {}
spacing['drift'] = 20.0 # m
spacing['pkg'] = 20.0 # m
# The EBS design for the salt repository requires some discussion. Because the
# thermal conductivity of crushed salt is several times less than for intact salt,
# the calculation radius for the homogeneous calculation is set at 4 m, somewhat
# farther than the 3.048 m radius if the backfill is converted, volumetrically,
# to a cylindrical geometry. Heat-generating waste packages would be placed into
# semi-cylindrical cavities milled in the alcove floor to improve heat transfer
# to the intact salt; thus, at least half of the waste package circumference is
# in close contact with intact salt. The reference uses intact salt properties
# from 4 m inward to the waste package, but with only 75% of the periphery
# available to transfer heat. The combination of half the package surface in
# contact with intact salt, and half with backfill, is represented by 75% in
# contact with intact salt.
# Because the dimensions are based on thickness, a backfill thickness of 3.5 m is
# used to approximate the amount of crushed salt around any waste canister. HLW
# packages will have a radii of ~0.7 m, while SNF packages will have radii of
# approximately 0.8-1.2 m
# canister
# EBS dimensions for spent nuclear fuel:
# Backfill thickness based on expected package radius ~ 0.41 m
ebs = {}
ebs['snf'] = {}
ebs['snf']['layer'] = ['backfill', 'overpack']
ebs['snf']['material'] = ['crushed salt', 'carbon steel']
ebs['snf']['dr'] = [3.590, 0.05] # calculation radius = 4.0 m
# SNF package diameter is 0.82 m
# SNF canister diameter is 0.72 m
# SNF overpack thickness is 0.05 m
# EBS dimensions for high level waste - no overpack
# Backfill thickness based on expected canister r = 0.305 m
ebs['hlw'] = {}
ebs['hlw']['layer'] = ['backfill', 'overpack']
ebs['hlw']['material'] = ['crushed salt', 'carbon steel']
ebs['hlw']['dr'] = [3.695, 0.0] # calculation radius = 4.0 m
###################################################################################