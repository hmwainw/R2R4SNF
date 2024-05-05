###############################################################################
"""
    Last edited on June 20, 2019

    @author: matz

    comments: Control MCNP iterations to determine the critical volume (cm^3)
              for a reflected sphere of rock, water, and heavy metal
              
"""
###############################################################################
import os
import numpy as np
from decimal import Decimal
from scipy.interpolate import interp1d
from scipy.optimize import root
import mcnp
import __init__
###############################################################################


class CriticalVolume(object):


    """
    Class to iterate on MCNP I/O to determine the critical volume (cm^3) 
    from the critical radius (cm) for a reflected spherical mixture of rock, 
    water, and heavy metal representing a deposition of fissile material in 
    the far field of a geologic repository.
    
    The iteration scheme is as follows:
    1.  The infinite multiplication factor is calculated using MCNP. If the
        sphere is subcritical, the k-effective will be less than one even
        with the reflective boundary condition and the calculation will be
        terminated.
    2.  If the k-effective is greater than one, the code begins iterations 
        for finite spheres with a 1-m thick reflector made of water-saturated
        host rock. For the first iteration, the radius is guessed using an 
        analytical solution based on diffusion theory, with nuclear data 
        evaluated at 0.625 eV (the boundary between epithermal and thermal
        energies).
    3.  After the first iteration, the radius is updated by guessing based
        on the previous best result for k-effective. A small perturbation
        is applied so that if the result is not improved, the outcome of a
        subsequent iteration will be slightly different, allowing more points
        to be generated near the radius where k-effective = 1.
    4.  The final value of the critical radius, which is used to calculate
        the critical volume, is interpolated based on the data produced in
        the previous iterations.
    
    Currently, the neutronics calculations are performed using MCNP, but in 
    the future it would be ideal to perform this calculation with KENO in
    order to make the nwpy package only depedent on the SCALE package.
    
    """
    
    
    k_trgt = 0.98 # target k-effective


    ###########################################################################
    # CRITICAL MASS ITERATIONS
    # First two methods are the main operators
    ###########################################################################
    @classmethod
    def calculate(cls, name, core, core_dens, refl, refl_dens, wdpath='./',
                  tol=5e-3, verbose=False, print_mcnp=False):
        """Determine the critical volume of a sphere of SNM
        
        Parameters
        ----------
        name: str
            Case name for file I/O
            
        core: dict
            Composition of deposition core material by mass
        
        core_dens: float
            Density of the material described by core (g/cc)
            
        refl: dict
            Composition of deposition reflector material by mass
        
        refl_dens: float
            Density of the material described by refl (g/cc)
            
        wdpath (optional): str
            Working directory for MCNP I/O
        
        tol (optional): float
            Solution tolerance for the k-effective calculation
            
        verbose (optional): Boolean
            Indicates whether to print results with each iteration
        
        print_mcnp (optional): Boolean
            Indicates whether to print MCNP output to screen during runtime

        Returns
        -------
        Critical volume of deposition (cm^3)
        
        """
    
        if(not os.path.exists(wdpath)):
            os.makedirs(wdpath)
        r_crit = cls.calculate_radius(name, core, core_dens, refl,
                                      refl_dens, wdpath)
        if(r_crit is not None):
            return((4/3.0)*np.pi*r_crit**3) # cm^3
        else: # criticality screened
            return
    
    
    @classmethod
    def calculate_radius(cls, name, core_mat, core_dens, refl_mat, refl_dens,
                         wdpath='./', max_iter=10, tol=1e-2, verbose=False,
                         print_mcnp=False):
        """Determine the critical radius of a sphere of SNM
        
        Parameters
        ----------
        name: str
            Case name for MCNP file I/O
            
        core_mat: dict
            Composition of deposition core material by mass
            
        core_dens: float
            Density of material described by core_mat (g/cc)
        
        refl_mat: dict
            Composition of reflector material by mass
        
        refl_dens: float
            Density of material described by refl_mat (g/cc)
            
        wdpath (optional): str
            Working directory for MCNP I/O
        
        max_iter (optional): int
            Maximum number of MCNP calculation iterations
         
        tol (optional): float
            Solution tolerance for the k-effective calculation
            
        verbose (optional): Boolean
            Indicates whether to print results with each iteration
        
        print_mcnp (optional): Boolean
            Indicates whether to print MCNP output to screen during runtime
            
        Returns
        -------
        Critical radius of sphere (cm)
        
        """
        
        cls.k = [] # list of keff values
        cls.sd = [] # list of standard deviations from mcnp
        cls.r = [] # list of sphere radii
        cls.iter = 0 # iteration counter
        if(not os.path.exists(wdpath)):
            os.makedirs(wdpath)
        # First, do k-infinite calculation
        k_inf, sd = cls.mcnp_iteration(name+'_inf', wdpath, core_mat,
                                       core_dens, refl_mat, refl_dens)
        if(float(k_inf) < cls.k_trgt): # criticality screened
            return
        r0 = cls.guess_radius0(core_mat, core_dens, refl_mat, refl_dens, k_inf)
        error = 1.0
        while(error > tol and cls.iter < max_iter or cls.iter < 3):
            k, sd = cls.mcnp_iteration(name, wdpath, core_mat, core_dens,
                                       refl_mat, refl_dens, r0)
            if(verbose):
                cls.print_iteration(r0, k)
            r0, error = cls.update_iteration(r0, k, k_inf, sd)
        cr = cls._interpolate(cls.r, cls.k)
        cls.criticalradius = cr
        return(cr)
    
    
    ###########################################################################
    # MCNP I/O
    ###########################################################################
    @classmethod
    def mcnp_iteration(cls, name, wdpath, c_mat, c_dens, r_mat,
                       r_dens, r0=None, verbose=False):
        """Write the MCNP input, run MCNP, and retrieve the output"""
        
        if(r0 == None): # if iteration zero, do k-infinity calculation
            inp = cls.write_mcnp_input(name, wdpath, c_mat, c_dens,
                                       r_mat, r_dens, 500.0, True)
        else: # do k-eff calculation with input
            inp = cls.write_mcnp_input(name, wdpath, c_mat, c_dens,
                                       r_mat, r_dens, r0)
        out = mcnp.MCNP.run_mcnp(inp, verbose)
        cls.cleanup(inp) # delete src and run tapes
        k, sd = cls.read_cm_output(out)
        return(k, sd)


    @classmethod
    def write_mcnp_input(cls, name, wdpath, c_comp, c_dens, r_comp,
                         r_dens, radius, kinf=False):
        """Write the input file according to the iteration status"""
        
        input_filename = name+'_'+str(cls.iter)
        path_to_infile = os.path.join(wdpath, input_filename)
        mcnp.MCNP.make_input(path_to_infile, c_comp, c_dens, r_comp,
                             r_dens, radius, kinf)
        if(os.path.exists(path_to_infile+'o')):
            os.remove(path_to_infile+'o') # remove old output file
        return(path_to_infile)


    @staticmethod
    def cleanup(path_to_inputfile):
        """Clean up the source and run tapes"""
        
        os.remove(path_to_inputfile+'s')
        os.remove(path_to_inputfile+'r')


    @staticmethod
    def read_cm_output(path_to_outfile):
        """Read MCNP output file, returning keff and standard deviation"""
        
        output = open(path_to_outfile).read()
        tag = output.find('keff = ')
        keff = output[(tag+7):(tag+14)]
        sd = output[(tag+55):(tag+62)]
        return([keff, sd])


    ###########################################################################
    # ITERATION UPDATE
    ###########################################################################
    @classmethod
    def print_iteration(cls, radius, keff):
        """Print the output of the most recent iteration"""
        
        line = 'iter '+str(cls.iter)
        line += ' r='+str(Decimal(radius).quantize(Decimal('0.0001')))
        line += ' k='+str(Decimal(keff).quantize(Decimal('0.0001')))
        print(line)


    @classmethod
    def update_iteration(cls, radius, keff, kinf, err):
        """Update internal parameters based on most recent MCNP output"""
    
        cls.r.append(radius)
        cls.k.append(float(keff))
        cls.sd.append(float(err))
        best_k, best_r = cls._find_closest_keff()
        r0 = cls.guess_radius(best_k, best_r, float(kinf))
        error = abs(cls.k_trgt - cls._find_closest_keff()[0]) # best k
        cls.iter += 1
        return(r0, error)
    

    @classmethod
    def _find_closest_keff(cls):
        """Look through keff results, pick the one closest to unity"""
        
        t = [abs(x-cls.k_trgt) for x in cls.k]
        sorted_list = sorted(enumerate(t), key=lambda x:x[1])
        best_k = cls.k[sorted_list[0][0]]
        best_r = cls.r[cls.k.index(best_k)]
        return(best_k, best_r)
    
    
    @classmethod
    def _interpolate(cls, r, k):
        """Use linear interpolation to determine the critical radius"""
        
        interp = interp1d(k, r, bounds_error=False, fill_value='extrapolate')
        return(interp(cls.k_trgt))
    
    
    ###########################################################################
    # CALCULATE ITERATION 0 RADIUS
    # Take a guess at the radius of the first sphere by applying diffusion
    # theory, using cross section data of the major actinide isotopes in
    # at 1 MeV.
    ###########################################################################
    @classmethod
    def guess_radius0(cls, c_comp, c_dens, r_comp, r_dens, kinf):
        """From diffusion theory, take a guess at the initial radius for the
        calculation based on the composition of the material and the thermal
        cross sections. First, calculate the unreflected critical radius; 
        then, using that as a starting point, calculate the reflected radius
        """
        
        # unreflected radius - use as starting point for transcendental soln
        R_unrfl = cls.unreflected_radius0(c_comp, c_dens, float(kinf))
        # neutronic data for core and reflector: abs, tr, c, f, nuf
        # use to get Lsq
        #core = cls.calc_macro_xs(c_comp, c_dens, data.xs, data.mm, data.nu)
        #refl = cls.calc_macro_xs(r_comp, r_dens, data.xs, data.mm, data.nu)
        #D_core = 1/3.0/core[1]
        #B = np.pi/R_unrfl
        #D_refl = 1/3.0/refl[1]
        #L_refl = np.sqrt(D_refl/refl[0])
        #R = cls.solve(R_unrfl, (B, D_refl, D_core, L_refl,))
        #return(R)
        return(R_unrfl)
        
    
    @staticmethod
    def solve(x0, args):
        """Solve the transcendental equation for reflected radius with an
        initial guess as input and neutronic data as required arguments"""
    
        # dummy function for the transcendental eqn
        def f(r, b, dr, dc, lr):
            return(dc*b/np.tan(b*r) - dc/r + dr*(r+lr)/(r*lr))
        # use the scipy.optimize.root method to get result
        res = root(f, x0, args=args)
        return(res.x[0])
    

    @classmethod
    def unreflected_radius0(cls, comp, dens, kinf):
        abs,tr,c,f,nuf=cls.calc_macro_xs(comp,dens,data.xs,data.mm,data.nu)
        D = 1/3.0/tr
        Lsq = D/abs
        R_ext = np.pi*np.sqrt(cls.k_trgt*Lsq/(kinf-cls.k_trgt))
        R = R_ext - 2*D
        return(R)
    

    @classmethod
    def calc_macro_xs(cls, composition, dens, xs_dict, a_dict, nu_dict):
        """Calculate the macroscopic absorption and transport cross
            sections for a mixture of actinides"""
        
        Sig_c = 0.0
        Sig_f = 0.0
        nuSig_f = 0.0
        Sig_tr = 0.0
        comp = __init__.renormalize(composition)
        for nuc in list(set(comp.keys()) & set(xs_dict.keys())):
            if('f' in xs_dict[nuc].keys()):
                Sig_f += comp[nuc]*xs_dict[nuc]['f']/a_dict[nuc]
                nuSig_f += Sig_f*nu_dict[nuc]
            Sig_c += comp[nuc]*xs_dict[nuc]['c']/a_dict[nuc]
            xs_tr = xs_dict[nuc]['t']-xs_dict[nuc]['s']*2.0/3/a_dict[nuc]
            Sig_tr += comp[nuc]*xs_tr/a_dict[nuc]
        Sig_abs = (Sig_c+Sig_f)*6.02e23*1e-24
        nuSig_f *= dens*6.02e23*1e-24
        Sig_tr *= dens*6.02e23*1e-24
        Sig_c *= dens*6.02e23*1e-24
        Sig_f *= dens*6.02e23*1e-24
        return(Sig_abs, Sig_tr, Sig_c, Sig_f, nuSig_f)

    
#    @classmethod
#    def calc_macroscopic_xs(cls, composition, dens, xs_dict, a_dict):
#        """Calculate the macroscopic absorption and transport cross
#        sections for a mixture of actinides"""
#    
#        Sig_abs = 0.0
#        Sig_tr = 0.0
#        comp = cls.renormalize(composition)
#        for nuc in list(set(comp.keys()) & set(xs_dict.keys())):
#            for rxn in [d for d in xs_dict[nuc].keys() if d in ['c', 'f']]:
#                Sig_abs += comp[nuc]*xs_dict[nuc][rxn]/a_dict[nuc]
#            xs_tr = xs_dict[nuc]['t']-xs_dict[nuc]['s']*2.0/3/a_dict[nuc]
#            Sig_tr += comp[nuc]*xs_tr/a_dict[nuc]
#        Sig_abs *= dens*6.02e23*1e-24
#        Sig_tr *= dens*6.02e23*1e-24
#        return(Sig_abs, Sig_tr)


    @staticmethod
    def renormalize(comp):
        """Renormalize a dictionary so that the values sum to one"""
        
        total = sum(comp.values())
        return(dict((k,v/total) for k,v in comp.iteritems()))
    
    
    ###########################################################################
    # CALCULATE ITERATION 1 RADIUS
    # The second iteration requires a guess for the radius based on the result
    # from the first calculation.
    ###########################################################################
    @classmethod
    def guess_radius(cls, k, r, kinf):
        """Calculate the radius after the first iteration; update the value
        for the diffusion area based on the new knowledge of k0 and r0 and
        calculate a new value for the radius
        
        k: float
            Previous best result for k-effective
            
        r: float
            Radius corresponding to best result for k-effective
        
        kinf: float
            k-infinity from the screening calculation
        
        """
        
        delta = 2.0 # extrapolation distance
        if(cls.iter > 0):
            perturbation = (0.25+np.random.random())/float(cls.iter+1)
        else:
            perturbation = 0.0
        Bgsq = (np.pi/(r+delta))**2
        Lsq = ((kinf/k)-1)/Bgsq
        return(np.pi*np.sqrt(Lsq/(kinf/cls.k_trgt-1))-delta-perturbation)


################################################################################







