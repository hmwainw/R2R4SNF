####################################################################################################
"""
	Last modified: November 30, 2017
	
	@author: Milos Atz <milos.atz@berkeley.edu>
	
	comments: Run file for repository footprint calculation
	
"""
####################################################################################################
# INPUT FOR HEAT GENERATION CODE
####################################################################################################
# WASTE TYPE
waste_type = ['uox', 'mox', 'coex', 'nuex', 'ecc', 'ecm']
#waste_type = 'uox'
####################################################################################################
# HOST ROCK TYPE
rock_type = ['granite', 'clay', 'salt']
#rock_type = 'granite'
####################################################################################################
# TIME
# Surface storage time before emplacement.
storage_time = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]	# y
####################################################################################################
import numpy as np
import os
from optimf import optimf
apc = {}
for r in range(0, len(rock_type)):
	rock = rock_type[r]
	# define np array with ncols = n_wastetypes and nrow = n_storagetpts
	apc[rock] = np.zeros((len(storage_time), len(waste_type)))
	for w in range(0, len(waste_type)):
		waste = waste_type[w]
		for t in range(0, len(storage_time)):
			# count backwards from largest storage time
			st = storage_time[t]
			# if the previous two cases are equal, the min has been reached, don't run anymore
			# just set all of the subsequent values to be equal
			if(t>2 and apc[rock][t-1,w]>0 and round(apc[rock][t-2,w],2)==round(apc[rock][t-1,w],2)):
				apc[rock][t:len(storage_time),w]=apc[rock][t-1,w]
				break
			else:
				tmp = optimf(rock, waste, st, verbose=True)
				#tmp = np.random.random()
				# optimf returns nothing if constraint cannot be met
				if(tmp):
					apc[rock][t,w] = tmp
	if(not os.path.exists('../output')):
		os.mkdir('../output/')
	dat = open('../output/'+rock+'.apc' ,'w')
	l = 'st,'+','.join(map(str,waste_type))+'\n'
	for t in range(0, len(apc[rock])):
		l = l+str(storage_time[t])+','+','.join(map(str,apc[rock][t]))+'\n'
	dat.writelines(l)
	dat.close()
####################################################################################################
