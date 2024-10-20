import os

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_Ref_PWR')
os.system('~/pflotran/src/pflotran/pflotran >log &')

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_SPWR')
os.system('~/pflotran/src/pflotran/pflotran >log &')

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_HTGR')
os.system('~/pflotran/src/pflotran/pflotran >log &')

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_HPR')
os.system('~/pflotran/src/pflotran/pflotran >log &')

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_HTGR_FCM')
os.system('~/pflotran/src/pflotran/pflotran >log &')

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_SFR')
os.system('~/pflotran/src/pflotran/pflotran >log &')

