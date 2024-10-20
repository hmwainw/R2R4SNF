###################################################################################
"""
    Last modified: August 13, 2019
    
    @author: Milos Atz <milos.atz@berkeley.edu
    
    comments: Class controlling the iterations to calculate repository footprint
    
"""
###################################################################################
import os
import datetime
import numpy as np
###################################################################################


class Iteration(object):


    """
    ITERATION LOGGING AND CONTROL
    The Iteration class contains methods that control and log the iterations in
    the calculation to determine the repository footprint. 
    """
    
    
    def __init__(self):
        self.array_idx = 0
        self.iter_idx = 0
        self.data = {}
        self.data['wps'] = []
        self.data['ds'] = []
        self.data['area'] = []
        self.data['temp'] = []


    def __repr__(self):
        """repr method for Iterator class"""
        
        return('Iterator object')
    
    def __str__(self):
        """str method for Iterator class"""
        
        # assume max possible iterations < 999; len is 3
        # max are is on the order of 100,000.000; len is 9
        
        p = 'iter '+str(self.iter_idx)+';'+' '*(3+1-len(str(self.iter_idx)))
        p += 'DS = '+"{0:.3f}".format(self.data['ds'][-1])+'  '
        p += 'WPS = '+"{0:.3f}".format(self.data['wps'][-1])+'  '
        area = "{0:.3f}".format(self.data['area'][-1])
        p += 'A = '+area+' '*(10+1-len(area)) # spacing after area
        p += 'WP Peak Temp = '+"{0:.3f}".format(self.data['temp'][-1])
        return(p)
    

    def reset(self):
        """Reset the internal data for a new calculation"""
        
        self.array_idx += 1
        self.iter_idx = 0
        self.data['wps'] = []
        self.data['ds'] = []
        self.data['area'] = []
        self.data['temp'] = []
    
    
    def update(self):
        """Update the iteration index"""
        
        self.iter_idx += 1
        
        
    def read(self, **kwargs):
        """Update internal data with each footprint calc iteration"""
                
        for key in kwargs.keys():
            try:
                self.data[key] = np.append(self.data[key], kwargs[key])
            except:
                raise


    def get_peak_temp(self, area):
        """Find the peak temperature given the minimized area"""
        
        idx = np.where(self.data['area']==area)[0][0]
        return(self.data['temp'][idx])


    def write_data(self, rep, res, sens, sens_limit=0.05, rw_temp=None):
        """Write the data from a series of footprint iterations 
        with a fixed N to a data file
        
        Parameters
        ----------
        self: iter instance
            Contains iteration history of array dimensions and peak temperature
        
        rep: Array instance
            Repository array instance containing data about the problem
            
        res: scipy.optimize.minimze result
            Result of footprint minimization constrained by peak temperature
        
        sens: float
            Sensitivity of the solution to increasing the array size to N+2 x N+2
        
        rw_temp: float
            Temperature (degrees C) of the rock wall at the minimized dimensions
        
        """
        
        if(self.array_idx==0):
            # open the file as "write"
            openmode = 'w'
            # write the header
            t =  'Repository Footprint Minimization Code v.0.1 \n'
            t += 'Subpackage in the NWPY Python Package \n'
            t += 'Written by Milos Atz \n'
            t += 'Nuclear Waste Management Group \n'
            t += 'University of California, Berkeley \n'
            # write the time
            t += ('Output file generated: '+
                  str(datetime.datetime.today())[:-7]+'\n')
            # write the run information
            t += ('#######################################'+
                  '######################################\n')
            t += 'Host rock, wasteform, and repository data \n'
            t += ('Limit for effect of adjacent sources: '+
                  str(round(100*sens_limit))+'% \n')
            t += 'Repository type: '+str(rep.name)+'\n'
            t += 'Rock thermal conductivity [W/m/K]: '+str(rep.rock['k'])+'\n'
            t += 'Rock thermal diffusivity [m^2/s]: '+str(rep.rock['a'])+'\n'
            t += 'Waste ID: '+str(rep.waste['id'])+'\n'
            try:
                t += 'Waste form: '+str(rep.waste['form'])+'\n'
            except:
                pass
            t += 'Package loading (units): '+str(rep.waste['n_wf'])+'\n'
            t += 'Waste package length [m]: '+str(rep.pkg['l'])+'\n'
            t += 'Waste package diameter [m]: '+str(rep.pkg['d'])+'\n'
            t += 'Surface storage time [y]: '+str(rep.st)+'\n'
            t += ('Repository ambient temperature [C]: '+
                  str(rep.ambient_temp)+'\n')
            t += ('Waste package temperature limit [C]: '+
                  str(rep.constraint['wp'])+'\n')
            t += ('Rock wall temperature limit [C]: '+
                  str(rep.constraint['rock'])+'\n')
        else:
            # open the outfile as "append"
            openmode = 'a'
            t=''
        # write the output data from the optimization cycles
        t += ('#######################################'+
              '######################################\n')
        t += ('Optimization run history and convergence: '+\
              str(rep.N)+'x'+str(rep.N)+' array'+'\n')
        t += ('iter'+'\t\t'+
              'drift spacing'+'\t'+
              'wp spacing'+'\t'+
              'footprint'+'\t'+
              'wp peak temp'+'\n')
        for i in range(0, self.iter_idx):
            t += (str(i)+'\t\t'+
                  "{0:.4f}".format(self.data['ds'][i])+'\t\t'+\
                  "{0:.4f}".format(self.data['wps'][i])+'\t\t'+\
                  "{0:.4f}".format(self.data['area'][i])+'\t'+\
                  "{0:.4f}".format(self.data['temp'][i])+'\n')
        t += str(res)+'\n'
        t += ('#######################################'+
              '######################################\n')
        t += 'Constraints and sensitivity\n'
        if(rw_temp is not None):
            t += ('Rock wall temperature at current dimensions [C]: '+
                  "{0:.4f}".format(rw_temp)+'\n')
        t += ('Relative impact of increasing array size to '+
              str(rep.N)+'x'+str(rep.N)+': '+"{0:.2f}".format(100*sens)+'%\n')
        # if the path to the dir doesn't exist, make the dir
        if(not os.path.isdir(rep.wdpath)):
            os.makedirs(rep.wdpath)
        outname = os.path.join(rep.wdpath, self._make_outfile_name(rep))
        outfile=open(outname, openmode)
        outfile.write(t)
        outfile.close()


    @staticmethod
    def _make_outfile_name(rep_inst):
        """Define output file name using Repository inst attributes"""
        
        outname = rep_inst.waste['id']+'_'
        if('form' in rep_inst.waste.keys()):
            if(rep_inst.waste['form'] != rep_inst.waste['id']):
                outname += rep_inst.waste['form']+'_'
        outname += 'st='+str(rep_inst.st)+'_'
        outname += 'nwf='+str(rep_inst.waste['n_wf'])+'.out'
        if(rep_inst.waste['id']=='hlw' and
           'loading_fraction' in rep_inst.waste.keys()):
            outname += '_lf='+str(rep_inst.waste['loading_fraction'])
        return(outname)


#################################################################################
