###############################################################################
"""
    Last edited on August 31, 2018

    @author: matz
    
    comments: Script to benchmark all FCES Evaluation Group Analysis Examples

"""
###############################################################################
from nwpy.fuelcycle import fuelcycle
###############################################################################
#no_recycle = ['01','02','03','04','05','06','07','08']
#lim_recycle = ['09','10','11','12','13','14','15','16','17','18']
cont_recycle = ['19','20','21','22','23','24','25','26','27','28','29',
                '30','31','32','33','34','35','36','37','38','39','40']
###############################################################################
# Make this a pytest function so that I can call with pytest from main dir
# For example, have benchmark return a tuple with the three values
# Can also give benchmark a "print" option to print, like it already does
#fces = no_recycle+lim_recycle+cont_recycle
for num in cont_recycle:#fces:
    name = 'eg'+num
    print(name.upper())
    try:
        fc = fuelcycle.FuelCycle(name)
    except:
        print('fuel cycle '+name+' initialization failed')
        continue
    else:
        fc.benchmark()


###############################################################################
