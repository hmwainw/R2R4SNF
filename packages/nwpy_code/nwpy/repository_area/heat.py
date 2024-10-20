###################################################################################
"""
	Last modified: August 27, 2018
	
	@author: Milos Atz <milos.atz@berkeley.edu
	
	comments: Class methods for heat transfer calculations for repository array
	
"""
###################################################################################
import numpy as np
import math as math
import warnings
from scipy.integrate import quad
###################################################################################


class HeatTransfer(object):
    
    
    """
    This class contains the methods necessary to evaluate the temperature 
    constraint for an array of heat-generating waste packages in a close-contact
    geological repository. More information on the models behind this code can be
    found in reports by Sandia National Laboratories, namely:
    
    E. Hardin et. al. "Generic Repository Design Concepts and Thermal Analysis"
    SAND2011-6202, Sandia National Laboratories, 2011
    
    """


    def __init__(self, st=0.0, tol=1.0e-4):
        self.st = st
        self.tol=tol


    #------------------------------------------------------------------------------
    # PEAK ROCK TEMPERATURE
    # This function returns the temperature at the interface between the drift
    # wall and the EBS.
    #------------------------------------------------------------------------------
    def find_peak_rock_temp(self, rep_inst):
        """Evaluate the peak temperature of the rock
            
        Parameters
        ----------
        rep_inst: Repository instance 
            contains rock and waste data
            
        Returns
        -------
        Peak temperature of the rock (degrees C)
            
        """
                
        T_wall = np.array([0.0])
        t0 = 0.01 # first time point (years after disposal)
        T_wall = np.append(T_wall, self._calc_wall_temp(t0, rep_inst))
        step = 0
        while(T_wall[step+1]-T_wall[step] >= 0):
            step+=1
            t=10**(np.log10(t0)+step*0.1)
            T_wall = np.append(T_wall, self._calc_wall_temp(t, rep_inst))
        return(max(T_wall))
    
    
    #------------------------------------------------------------------------------
    # PEAK WASTE PACKAGE TEMPERATURE
    # Given the temperature history at the drift wall, the temperature through the
    # EBS layers can be calculated. For every point in the drift wall temperature
    # history, this function returns the temperature at the surface of the central
    # waste package.
    #------------------------------------------------------------------------------
    def find_peak_wp_temp(self, r, allsources=True):
        """Calculate the maximum temperature at the surface of the waste package
        do to its decay heat and the decay heat from adjacent heat sources.
        
        Parameters
        ----------
        
        r: Repository instance 
            contains rock and waste data
        
        allsources(optional): Boolean 
            Indicates whether to calculate the peak waste package temperature 
            based on all heat in array sources or just the central package
        
        Returns
        -------
        Peak temperature of the waste package surface (degrees C)
        
        """
            
        T_wps = np.array([0.0])
        # evaluate first time point
        t0 = 0.01 # first time point (years after disposal)
        T_wall = self._calc_wall_temp(t0, r, allsources)
        T_wps = np.append(T_wps, self._calc_wp_temp(t0, r, T_wall))
        step = 0
        while(T_wps[step+1]-T_wps[step] >= 0 or step < 20):
            # from the first point evaluated, temperature may initially
            # decrease because the heat from adjacent packages will take
            # some time to reach the waste package. When the waste package
            # is initially hot enough (like with low surface storage time, it
            # will begin to dissipate heat before the temperature contribution
            # from adjacent packages is accounted for. Therefore, in order to
            # definitely capture the peak temperature, a minimum number of points
            # is required in order to capture the peak temperature.
            step+=1
            t=10**(np.log10(t0)+step*0.1)
            T_wall = self._calc_wall_temp(t, r, allsources)
            T_wps = np.append(T_wps, self._calc_wp_temp(t, r, T_wall))
        return(max(T_wps))


    #------------------------------------------------------------------------------
    # PEAK SINGLE WASTE PACKAGE TEMPERATURE
    # If the disposal of a single package violates the repository temperature
    # constraints, there is no feasible solution. This function evaluates whether
    # the disposal of a single canister violates that constraint and returns the
    # corresponding Boolean.
    #------------------------------------------------------------------------------
    def single_package(self, rep_inst):
        """Check if disposal of a single package violates temperature constraints
        
        Parameters
        ----------
        rep_inst: Repository instance 
            contains rock and waste data
        
        Returns
        -------
        Boolean indicating if disposal of a single package violates surface 
        temperature limit
        
        """
            
        T_peak = self.find_peak_wp_temp(rep_inst, allsources=False)
        if(T_peak > rep_inst.constraint['wp']):
            return(False)
        else:
            return(True)


    #------------------------------------------------------------------------------
    # SENSITIVITY TO ARRAY SIZE
    # When the spacing between the drifts and packages is decreased, the peak
    # waste package surface temperature may be sensitive to the implied canisters
    # outside the simulated array due to the decrease in array size. This function
    # calculates the change in waste package peak temperature at current
    # repository dimensions if an additional row and column of waste packages are
    # added.
    #------------------------------------------------------------------------------
    def calc_array_size_sensitivity(self, rep_inst, T_peak_old):
        """Evaluate the sensitivity of peak WP temp to increased array size
        
        Parameters
        ----------
        rep_inst: Repository instance
            contains rock and waste data
        
        T_peak_old: float
            Peak waste package surface temperature at the initial array size
            
        Returns
        -------
        Float indicating the relative error in the peak temperature incurred by 
        increasing the repository array size.
        
        """
            
        # Add 2 extra drifts, as well as two  extra packages to every drift
        # i.e. 7x7 -> 9x9
        rep_inst.N += 2
        # Calculate new WP peak temp
        T_peak_new = self.find_peak_wp_temp(rep_inst)
        # undo this - otherwise would be too hidden
        rep_inst.N += -2
        # Evaluate relative to "old" for conservatism ("old" is lower temp)
        return(abs(T_peak_new-T_peak_old)/T_peak_old)


    #------------------------------------------------------------------------------
    # CALCULATION OF DRIFT WALL TEMPERATURE
    # All canisters in the repository array will contribute heat to the evaluation
    # point, which is taken as the drift wall in the center of the repository. In
    # this function, the contribution of heat from each canister is accounted for
    # to determine the temperature at the drift wall for some time after disposal.
    # Some symmetry is employed to decrease calculation time.
    #------------------------------------------------------------------------------
    def _calc_wall_temp(self, time, r, allsources=True):
        """Calculate temperature increase at the evaluation point from adjacent 
        heat sources
        
        Parameters
        ----------
        time: Time after emplacement in repository (years)
        
        r: Repository instance (contains rock and waste data)
        
        allsources (optional): Boolean 
            Indicating whether to return the wall evaluated based on all heat 
            sources in array or just the central package
        
        Returns
        -------
        Temperature history at the drift wall
        
        """
            
        T_out = {}
        # Finite line contribution
        T_out['fl'] = self._finiteline(time, r)
        if(not allsources):
            return(T_out['fl']+r.ambient_temp)
        else: # point contributions from other canisters in array
            T_out['pt']=0.0
            for drft in range(0,int((r.N-1)/2.0)+1):
                y = drft*r.spacing['drift']
                if(drft==0): # central drift treated differently
                    for can in range(1, int((r.N+1)/2.0)):
                        x = can*r.spacing['pkg']
                        d = np.sqrt(x**2+y**2) # distance between x and y
                        T_out['pt']+=2*self._point(time,r,d)
                else:
                    for can in range(0, int((r.N+1)/2.0)):
                        x = can*r.spacing['pkg']
                        d = np.sqrt(x**2+y**2) # distance between x and y
                        if(can==0): # first can treated differently
                            T_out['pt']+=2*self._point(time,r,d)
                        else:
                            T_out['pt']+=4*self._point(time,r,d)
            return(T_out['pt']+T_out['fl']+r.ambient_temp)


    #------------------------------------------------------------------------------
    # CALCULATION OF WASTE PACKAGE SURFACE TEMPERATURE
    # Given the temperature at the drift wall, the temperature through the EBS
    # layers can be calculated by steady-state cylindrical conduction (the SS
    # assumption is appropriate due to the relatively low thermal mass of the EBS
    # materials; it implies that heat transfer through the EBS is relatively fast
    # relative to heat transfer through the surrounding rock). This function
    # returns the temperature at the surface of the central waste package.
    #------------------------------------------------------------------------------
    def _calc_wp_temp(self, time, r, T_out):
        """Calculate the temperature across concentric EBS layers via conduction
        
        Parameters
        ----------
        time: float
            Time after emplacement in repository (years)

        r: Repository instance
            contains rock and waste data
        
        T_wall: float
            Temperature at the evaluation point from outside model (deg C)
        
        Returns
        -------
        Temperature at the waste package surface
        
        """
            
        T_ebs = {}
        qL = r.decay_heat(time+self.st)/r.pkg['l'] # length-averaged heat
        # Iterate over EBS layers to calculate inner radial temperature; the last
        # entry in the ebs dict is "wp" for waste package. The while loop should
        # quit when it sees that, at which point the last temp calculated will be
        # that at the surface of the waste package.
        layer_indx= 0
        r_out = r.ebs['r_drift']
        for l in r.ebs['layer']:
            if(l != 'overpack'):
                dr = r.ebs['dr'][layer_indx]
                r_in = r_out-dr
                k = r.ebs['k'][layer_indx]
                T_ebs[l] = T_out+qL*(np.log(r_out/r_in)/k)/2.0/np.pi
                # update the layer and "outer" temp for the next calculation
                T_out = T_ebs[l]
                layer_indx = layer_indx+1
                r_out = r_in
        # WP surface temp equals temp at inside radius of the last EBS layer
        last_layer = r.ebs['layer'][layer_indx-1] # just step back one
        return(T_ebs[last_layer])
    
    
    #------------------------------------------------------------------------------
    # HEAT CONTRIBUTION FROM FINITE LINE SOURCE
    # The central canister in the array is represented as a finite line source due
    # to its proximity to the evaluation point at the center of the array. This
    # function calculates the change in temperature over time at a fixed distance
    # (the drift wall) away from a time-dependent finite line source.
    #------------------------------------------------------------------------------
    def _finiteline(self, trgt, r):
        """Calculate temperature increase at distant point due to time variant
        finite line source centered at the origin
            
        Parameters
        ----------
        trgt: float
            Time after emplacement in repository (years)
        
        r: Repository instance 
            contains rock and waste data
            
        Returns
        -------
        Temperature increase due to finite line source
        
        """

        a = r.rock['a']*3600*24*365 # Convert time from seconds to years
        # integrand function
        def integrand(t):
            denom = 8.0*np.pi*r.rock['k']*(trgt-t)
            # (discarge time) = (storage time) + (relative time)
            hterm = r.decay_heat(t+self.st)/r.pkg['l']
            expterm = np.exp(-(r.ebs['r_drift']**2)/4.0/a/(trgt-t))
            erf1 = math.erf(0.5*(0.5*r.pkg['l'])/np.sqrt(a*(trgt-t)))
            erf2 = math.erf(0.5*(-0.5*r.pkg['l'])/np.sqrt(a*(trgt-t)))
            return(hterm*expterm*(erf1-erf2)/denom)
        # compute the integral
        return(self._integrate(integrand, trgt))


    #------------------------------------------------------------------------------
    # HEAT CONTRIBUTION FROM POINT SOURCE
    # The canisters adjacent to the central canister are represented as point
    # sources. This function calculates the change in temperature over time at
    # some distance away from a time-dependent point source.
    #------------------------------------------------------------------------------
    def _point(self, trgt, r, dist):
        """Calculate temperature increase at distant point due to a time variant 
        point source located at the origin
        
        Parameters
        ----------
        trgt: float
            Time after emplacement in repository (years)

        r: Repository instance 
            contains rock and waste data
        
        dist: float
            Center-to-center distance (m) between source and central canister
        
        Returns
        -------
        Temperature increase due to point source
        
        """
            
        a = r.rock['a']*3600*24*365 # Convert time from seconds to years
        total_dist = np.sqrt(r.ebs['r_drift']**2+dist**2)
        # integrand function
        def integrand(t,d):
            denom = 8.0*r.rock['k']*np.sqrt(a)*(np.pi**1.5)*((trgt-t)**1.5)
            expterm = np.exp(-(d**2)/4.0/a/(trgt-t))
            # (discarge time) = (storage time) + (relative time)
            hterm = r.decay_heat(t+self.st)
            return(hterm*expterm/denom)
        # compute the integral
        return(self._integrate(integrand, trgt, arguments=(total_dist,)))


    def _integrate(self, integrand, target, arguments=()):
        """Integrate a function from 0 to some value using 
        the scipy.interpolate.quad function"""
        
        abserr = 1.0
        counter = 1
        # NOTE: The number of integration subdivisions is automatically
        # increased if scipy.integrate.quad fails. To prevent warnings
        # every time, all warnings are ignored; this may be bad practice,
        # suggestions welcome.
        warnings.filterwarnings("ignore")
        while(abserr > self.tol): # Iterate to compute the integral
            x, abserr = quad(integrand, 0, target, args=arguments,
                             limit=50*counter)[0:2]
            counter=counter+1
            if(counter > 5):
                break
        # Reset warnings, check for convergence
        warnings.resetwarnings()
        if(abserr > self.tol):
            warnings.warn('Integral not converged after 250 subdivisions; +'\
                          'abserr = '+str(abserr),RuntimeWarning)
        return(x)
###################################################################################