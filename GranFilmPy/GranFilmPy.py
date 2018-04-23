DEBUG="False"

# --- Import libraries and files from the python interface

import os as os
import tempfile
import shutil
import matplotlib.pyplot as plt
import numpy as np

import GranFilm_parameters as gf_param
import GranFilm_results
import GranFilm_tools
from platform import system

operating_system = system()

# --- Set practical paths

Python_Files       = os.path.dirname(os.path.realpath(__file__)).replace('\\','/') # Path to the directory of the .py files constituting the python interface for GranFilm 
GranFilmPy_Package = os.path.split(Python_Files)[0].replace('\\','/') # Path to the GranFilmPy package directory
GRANFILM_ROOT      = os.path.join(GranFilmPy_Package,'bin','GranFilm_'+operating_system).replace('\\','/') # Path to the Fortran executable that will be called during the calculation
SOPRA_ROOT         = os.path.join(GranFilmPy_Package,'SOPRA_DataBase/').replace('\\','/') # Path to the SOPRA database which contains the dielectric functions of numerous materials

  
class GranFilm:

    """
    
    Main class of the GranFilmPy package, a python interface to GranFilm, a software for the
    simulation of the optical properties of supported thin granular films.
    An instance of this class represents a simulation.
    
    ----- Attributes -----
    
    param -> A dictionary inherited class which stores the main parameters of the simulation
             and creates from these an input file for the Fortran executable.
    
    ----- Methods -----
    
     __init__(self,GRANFILM_ROOT=GRANFILM_ROOT,SOPRA_ROOT=SOPRA_ROOT,radius=8.0,truncation_ratio=0.0,radius_ratios=[1.0],
              media=['air','mgo','ag','mgo'],theta=45.0,phi=0.0,energy_range=[1.5,5],pol='p',arrangement="Lattice",
              lattice_type="Square",island_island="None",lattice_constant=20.0,multip_pos_rat=0.0,number_en_points=300,
              multip_order=30,points_file="None",energy=[2.2],area_ratio_pot=2.0,number_pot_points=300,param_file=None):
              
    Initializes a new GranFilm simulation.
    
    Arguments :
        
        GRANFILM_ROOT (string)       : The path to the GranFilm executable (By default {}).
        
        SOPRA_ROOT (string)          : The path to the SOPRA database which contains the dielectric functions of many materials
                                       (by default {}).
        
        radius (float or float list) : The radius of the particles in nanometers. A float for a sphere or a list = 
                                       [radius parallel to the interface, radius perpendicular to the interface]
                                       for a spheroidal particle (default 8.0 nm sphere).

        truncation_ratio (float)     : (distance center of the particle - substrate) / radius
                                       Supported values: truncation_ratio in [0,1].
                                       0 if the center of the particle is at the interface
                                       1 if the particle is sitting at the top of the substrate
                                       (default 0.0).
                                              
        radius_ratios (float list)   : For a coated particle we can define for each layer the radius ratio as
                                       [radius of the layer] / [radius of the entire particle]. The radius ratios
                                       is then the list of the radius ratio from the outer shell to the core.
                                       For a non-coated particle radius_ratios [1.0] (default [1.0]).
    
        media (string list)          : A list of strings filled with the materials used for the simulation. media = ['m1','m2','m3','m4', ...]
                                       where m1 is the material from which the incident light comes, m2 is the substrate, m3 and m4
                                       correspond to the part of the outer shell which is above and below the substrate, and so on for
                                       the inner shells (default ['air','mgo','ag','mgo'] i.e a half sphere of silver on a mgo substrate).

        theta (float)                : The angle of incidence of the light source in degrees (default 45 degrees).
        
        energy_range (float list)    : A list with the minimum / maximum of the energy range (in eV) on which GranFilm will calculate
                                       the optical properties of the system. energy_range = [energy_min,energy_max] (default [1.5,5.0] eV).

        lattice_type (string)        : The geometry of the 2D lattice -> 'Square' or 'Hexagonal' (default 'Square').

        lattice_constant (float)     : The distance between two particles of the lattice in nanometers (defaut 20.0 nm).
        
        island_island (string)       : Indicates if the interaction bewtween the particles of the lattice has to be taken into account during the simulation.
                                       "None" if it is not the case, "Dipole" or "Quadrupole" will make each particle interacts with the others as
                                       a dipole / quadrupole (default "None").
        
        points_file (string)         : The path to the text file with all the points where GranFilm will calculate the electric potential.
                                       points_file = "None" if you don't want the potential to be, "classic" for the potential in a rectangle
                                       around the particle, "surface" if you want an estimation of the error on the potential (averaged
                                       on all the interfaces) and '.../path/to/my_points_file' if you want the potential at specific points
                                       (default 'None').
        
        area_ratio_pot (float)       : The ratio of the size of the rectangle in which the potential will be calculated
                                       (for points_file = 'classic') by the radius of the particles (default 2.0).
        
        energy (float list)          : The energies at which the potential is calculated in eV (default [2.2] eV).
    """

# --------------- Instance initialisation --------------- #


    def __init__(self,GRANFILM_ROOT=GRANFILM_ROOT,SOPRA_ROOT=SOPRA_ROOT,radius=8.0,truncation_ratio=0.0,radius_ratios=[1.0],
                 media=['air','mgo','ag','mgo'],theta=45.0,phi=0.0,energy_range=[1.5,5],pol='p',arrangement="Lattice",lattice_type="Square",
                 island_island="None",lattice_constant=20.0,multip_pos_rat=0.0,number_en_points=300,multip_order=30,points_file="None",
                 energy=[2.2],area_ratio_pot=2.0,number_pot_points=300,param_file=None):
        
        self._dict = {}
        
        # Init the dictionaries of the parameters

        param_dict = dict(GRANFILM_ROOT=GRANFILM_ROOT,SOPRA_ROOT=SOPRA_ROOT,radius=radius,truncation_ratio=truncation_ratio,radius_ratios=radius_ratios, media=media,
                          theta=theta,phi=phi, energy_range=energy_range, pol=pol,arrangement=arrangement,lattice_type=lattice_type, island_island=island_island,
                          lattice_constant=lattice_constant,multip_pos_rat=multip_pos_rat,number_en_points=number_en_points, multip_order=multip_order)

        potential_dict = dict(points_file=points_file,energy=energy,area_ratio_pot=area_ratio_pot,number_pot_points=number_pot_points)
        
        self.param = gf_param.GranFilm_Parameters(param_dict,potential_dict,param_file)
        

# --------------- Call method to run the simulation --------------- #

    def __call__(self):

        # ----- Write the points_file if needed (and set the fortran_radius) ----- #
        
        self.param.write_points_file()
        
        # ----- Write the parameters file ----- #
        
        tmpdir = tempfile.mkdtemp(suffix='_GranFilm')
        if (DEBUG=="True"): print ("Created directory %s ") % tmpdir
        
        # Define filenames
        input_file  = "Parameters.sif"
        output_file = "Results.dat"
        
        # Write the parameters file
        assert os.path.exists( tmpdir ), "Directory %s does not exist" % tmpdir   # is this needed
        filename = os.path.join(tmpdir,input_file)
        self.param.File( filename, "w",format="Fortran") 


        # ----- Run the simulations ----- #
        
        cmd = '%s -py -p %s -o %s' % (self.param["GRANFILM_ROOT"], input_file, output_file)
        if (DEBUG=="True"): print (cmd)
        out,err,retcode,time = GranFilm_tools.run_command( cmd, tmpdir )

        # Check for errors
        if (retcode != 0):
            print ('******STDOUT:******')
            print (out)
            print ('******STDERR:******')
            print (err)
            print ('*******************')
            print ('THERE WAS AN ERROR RUNNING GRANFILM')
            shutil.rmtree(tmpdir)
            return
               
        # ----- Read the simulation results ----- #
        
        data          = GranFilm_tools.get_results_from_file(  os.path.join(tmpdir,output_file)  )
        data_RefTrans = GranFilm_tools.get_results_from_file(  os.path.join(tmpdir,output_file + "_RefTrans")  )
        data_Epsilon  = GranFilm_tools.get_results_from_file(  os.path.join(tmpdir,output_file + "_Epsilon")  )

        data_Potential = None
        if self.param.Potential["points_file"] != 'None': # Add enventual data for the potential
            data_Potential = GranFilm_tools.get_results_from_file( os.path.join(tmpdir,"Results_potential.dat")  )
        
        # Clean up the un-needed directory structure
        shutil.rmtree(tmpdir)   
                
        # Populates the GranFilm instance with the results of the simulation  
        GranFilm_results.GranFilm_Results(self, data, data_RefTrans, data_Epsilon, data_Potential)
        
        # Print some final comments
        if (DEBUG=="True"): print ("GranFilm instance populated with simulation results!")

        return
    

# --------------- Method for the calculation of reflectivity --------------- #

    
    def Reflectivity(self,pol,var='R',angle='Default'):
        
        medium_1      = self.param["media"][0]
        medium_2      = self.param["media"][1]
        
        if angle == 'Default':
            theta = self.param["theta"]*np.pi/180.0
        else:
            theta = angle*np.pi/180.0
                
        y = None
        
        ### --- Calculation for each case
        
        if var == 'r': # Fresnel coefficients
            
            # Useful quantities
            
            n1,n2         = np.sqrt(self.Epsilon[medium_1]),np.sqrt(self.Epsilon[medium_2])                
            s1            = np.sin(theta)
            c1            = np.cos(theta) # cos(theta)
            c2            = np.cos(np.arcsin(s1*n1.real/n2.real)) # cos(theta_2)
            gamma         = self.Susceptibilities('gamma')*10
            beta          = self.Susceptibilities('beta')*10
            wavelength    = self.wavelength
            t_gamma       = 1j*gamma*(2*np.pi/wavelength) # = i(w/c)*gamma
            K1            = (n2*c1+n1*c2)*(1-0.25*(2*np.pi/wavelength)**2*gamma*beta*(n1*s1)**2)
            K2            = (n2*c1-n1*c2)*(1-0.25*(2*np.pi/wavelength)**2*gamma*beta*(n1*s1)**2)    
            
            if pol == 'p':
                y = (K2-1j*(2*np.pi/wavelength)*(gamma*c1*c2-n1**3*n2*beta*s1**2))/(K1-1j*(2*np.pi/wavelength)*(gamma*c1*c2+n1**3*n2*beta*s1**2))                
            elif pol == 's':
                y = (n1*c1-n2*c2+t_gamma)/(n1*c1+n2*c2-t_gamma)
                
        elif var == 'rf': # Fresnel coefficients for the flat system (no islands)
            
            # Useful quantities
            
            n1,n2         = np.sqrt(self.Epsilon[medium_1]),np.sqrt(self.Epsilon[medium_2])
            s1            = np.sin(theta)
            c1            = np.cos(theta) # cos(theta)
            c2            = np.cos(np.arcsin(s1*n1.real/n2.real)) # cos(theta_2)
            
            if pol == 'p':                
                y = (n2*c1-n1*c2)/(n2*c1+n1*c2)              
            elif pol == 's':                
                y = (n1*c1-n2*c2)/(n1*c1+n2*c2)  
                
        elif var == 'R': # Intensity coefficients
            
            if pol == 'p' or pol == 's':
                y = np.abs(self.Reflectivity(pol,'r',angle))**2                
            elif pol == 'ps':                
                y = (np.abs(self.Reflectivity('p','r',angle))**2 + np.abs(self.Reflectivity('s','r',angle))**2)/2.0
                
        elif var == 'Rf': # Intensity coefficients for the flat system (no islands)
            
            if pol == 'p' or pol == 's':                
                y = np.abs(self.Reflectivity(pol,'rf',angle))**2                
            elif pol == 'ps':                
                y = (np.abs(self.Reflectivity('p','rf',angle))**2 + np.abs(self.Reflectivity('s','rf',angle))**2)/2.0
                
        elif var == 'dR_R' and medium_1 != medium_2: # Diffential reflectivity (SDRS measurement)
            
            assert(pol in ['p','s','ps']), "ERROR: pol = {} is not a supported parameter.".format(pol)
            y = (np.abs(self.Reflectivity(pol,'R',angle)) - np.abs(self.Reflectivity(pol,'Rf',angle)))/np.abs(self.Reflectivity(pol,'Rf',angle))
            
        assert(y is not None), "ERROR: pol = {}, var = {} or angle = {} are not supported parameters.".format(pol,var,angle)
        
        return y
    
    
# --------------- Transmittance --------------- #

    
    def Transmittance(self,pol,var='T',angle='Default'):
        
        medium_1      = self.param["media"][0]
        medium_2      = self.param["media"][1]
        
        if angle == 'Default':
            theta = self.param["theta"]*np.pi/180.0
        else:
            theta = angle*np.pi/180.0
                
        y = None
        
        ### --- Calculation for each case
        
        if var == 't': # Fresnel coefficients
            
            # Useful quantities
            
            n1,n2         = np.sqrt(self.Epsilon[medium_1]),np.sqrt(self.Epsilon[medium_2])                
            s1            = np.sin(theta)
            c1            = np.cos(theta) # cos(theta)
            c2            = np.cos(np.arcsin(s1*n1.real/n2.real)) # cos(theta_2)
            gamma         = self.Susceptibilities('gamma')*10
            beta          = self.Susceptibilities('beta')*10
            wavelength    = self.wavelength
            t_gamma       = 1j*gamma*(2*np.pi/wavelength) # = i(w/c)*gamma
            K1            = (n2*c1+n1*c2)*(1-0.25*(2*np.pi/wavelength)**2*gamma*beta*(n1*s1)**2)    
            
            if pol == 'p':
                y = 2*n1*c1*(1+0.25*(2*np.pi/wavelength)**2*gamma*beta*(n1*s1)**2)/(K1-1j*(2*np.pi/wavelength)*(gamma*c1*c2+n1**3*n2*beta*s1**2))
                #N = 2*n1*c1*(1+0.25*(2*np.pi/wavelength))
            elif pol == 's':
                y = 2*n1*c1/(n1*c1+n2*c2-t_gamma)
                
        elif var == 'tf': # Fresnel coefficients for the flat system (no islands)
            
            # Useful quantities
            
            n1,n2         = np.sqrt(self.Epsilon[medium_1]),np.sqrt(self.Epsilon[medium_2])
            s1            = np.sin(theta)
            c1            = np.cos(theta) # cos(theta)
            c2            = np.cos(np.arcsin(s1*n1.real/n2.real)) # cos(theta_2)
            
            if pol == 'p':                
                y = (2*n1*c1)/(n2*c1+n1*c2)              
            elif pol == 's':                
                y = (2*n1*c1)/(n1*c1+n2*c2)  
                
        elif var == 'T': # Intensity coefficients
        
            # Useful quantities
            
            n1,n2         = np.sqrt(self.Epsilon[medium_1]),np.sqrt(self.Epsilon[medium_2])                
            c1            = np.cos(theta) # cos(theta)
            c2            = np.cos(np.arcsin(s1*n1.real/n2.real)) # cos(theta_2)
        
            if pol == 'p' or pol == 's':
                y = np.abs(n2/n1)*(c2/c1)*np.abs(self.Transmittance(pol,'t',angle))**2                
            elif pol == 'ps':                
                y = np.abs(n2/n1)*(c2/c1)*(np.abs(self.Transmittance('p','t',angle))**2 + np.abs(self.Transmittance('s','t',angle))**2)/2.0
                
        elif var == 'Tf': # Intensity coefficients for the flat system (no islands)
            
            if pol == 'p' or pol == 's':                
                y = np.abs(n2/n1)*(c2/c1)*np.abs(self.Transmittance(pol,'tf',angle))**2                
            elif pol == 'ps':                
                y = np.abs(n2/n1)*(c2/c1)*(np.abs(self.Transmittance('p','tf',angle))**2 + np.abs(self.Transmittance('s','tf',angle))**2)/2.0
                
        elif var == 'dT_T' and medium_1 != medium_2: # Diffential transmittance (SDRS measurement)
            
            assert(pol in ['p','s','ps']), "ERROR: pol = {} is not a supported parameter.".format(pol)
            y = (np.abs(self.Transmittance(pol,'T',angle)) - np.abs(self.Transmittance(pol,'Tf',angle)))/np.abs(self.Transmittance(pol,'Tf',angle))
            
        assert(y is not None), "ERROR: pol = {}, var = {} or angle = {} are not supported parameters.".format(pol,var,angle)
        
        return y