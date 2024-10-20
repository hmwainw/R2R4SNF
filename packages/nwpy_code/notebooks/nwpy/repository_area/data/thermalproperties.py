###################################################################################
"""
    Last edited on June 23, 2018
	
    @author: matz
	
    comments: Thermal properties of host rock and EBS materials for the 
              calculation of heat in a repository, as per the SNL/LLNL 
              2011 report.
	
"""
###################################################################################
# THERMAL CONDUCTIVITY (W/m-K)
k = {}
k['granite'] = 2.5 # 100 C
k['clay'] = 1.75 # 100 C
k['salt'] = 3.2 # 200 C
k['bentonite'] = 0.6 # 100 C, dry
k['steel'] = 46.0
k['carbon steel'] = 53.0
k['crushed salt'] = 0.75*k['salt'] # 200 C; see salt repository data for details
k['copper'] = 366.9
###################################################################################
# THERMAL DIFFUSIVITY (m^2/s)
a = {}
a['granite'] = 1.13e-6 # 100 C
a['clay'] = 6.45e-7 # 100 C
a['salt'] = 1.60e-6 # 100 C
###################################################################################

