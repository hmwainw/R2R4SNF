import os

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_Ref_PWR')
os.system('~/pflotran/src/pflotran/pflotran >log &')
    
os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_Micro_PWR')
os.system('~/pflotran/src/pflotran/pflotran >log &')

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_HTGR_80')
os.system('~/pflotran/src/pflotran/pflotran >log &')

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_eVinci')
os.system('~/pflotran/src/pflotran/pflotran >log &')

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_USNC')
os.system('~/pflotran/src/pflotran/pflotran >log &')

os.chdir('/home/hmw/Documents/NuclearWaste/chloe_PFLOTRAN/pflotran_SFR')
os.system('~/pflotran/src/pflotran/pflotran >log &')



