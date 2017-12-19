DEBUG="False"

import numpy as np
import os as os
import tempfile
import shutil
import matplotlib.pyplot as plt

import GranFilm_parameters as gf_param
import GranFilm_results
import GranFilm_tools
from platform import system

operating_system = system()

# --- Set essential paths

python_interface = os.path.dirname(os.path.realpath(__file__)).replace('\\','/') # Path to the directory of the .py files constituting the python interface for GranFilm 
gf_package = os.path.split(python_interface)[0].replace('\\','/') # Path to the GranFilm_Package directory
GRANFILM_ROOT= os.path.join(gf_package,'bin','GranFilm_'+operating_system).replace('\\','/') # Path to the GranFilm executable
SOPRA_ROOT = os.path.join(gf_package,'SOPRA_DataBase/').replace('\\','/') # Path to the SOPRA_DataBase which contains dielectric functions for numerous materials

  
class GranFilm:

    """
    :synopsis: Function which initializes a simulation with the configuration given by the user. |br| \
    **GRANFILM_ROOT** and **SOPRA_ROOT** are needed arguments. |br| \
    The other arguments are the parameters of the simulation, they are optionnal and\
    set with default values if not defined as arguments of the function.

    
                                
    :param GRANFILM_ROOT: the path to the GranFilm executable
    :type GRANFILM_ROOT: string

    :param SOPRA_ROOT: the path to the SOPRA database with all the .nk files
    :type SOPRA_ROOT: string
    
    :param radius: the radius of the particles. A float for a sphere or a list = [ radius parallel to the interface , radius perpendicular \
    to the interface ] for a spheroidal particle
    
    :type radius: float or float list

    :param truncation_ratio: **[distance center of the particle - substrate] / radius** |br| \
    ``Supported values: truncation_ratio in [0,1]`` |br| |br| \
    ``Example:`` |br| |br| \
    |image_truncation_ratio| |br| |br|    
    
    :type truncation_ratio: float
    
    :param radius_ratios: For a coated particle we can define for each layer the radius ratio as **[radius of the layer] / [radius of the particle]** |br| \
    The radius ratios is then the list of the radius ratio from the outer shell to the core. |br| \
    For a non-coated particle the radius ratios is 1.0 |br| |br| \
    ``Example:`` |br| |br| \
    |image_radius_ratios| |br| |br|
    
    :type radius_ratios: float or float list

    :param media: A string with the materials used for the simulation, separated by a coma -> media = "m1,m2,m3,m4, ..." |br| \
    m1 is the material from which the incident light comes, m2 is the substrate, m3 and m4 correspond to the part of the \
    outer shell which is above and below the substrate, and so on for the inner shells. |br| |br| \
    ``Note: The SOPRA database must contain the dielectric functions of the`` ``given materials.`` |br| |br| \
    ``Example:`` |br| |br| \
    |image_media| |br| |br|

    :type media: string

    :param theta: the angle of incidence
    :type theta: float

    :param energy_range: a list with the minimum / maximum of the energy range on which GranFilm will provide data -> energy_range = [energy_min,energy_max] |br| \
    ``Note: make sure that the dielectric functions of the materials used are`` ``defined over this energy range``
    :type energy_range: float list

    :param lattice_type: the type of the 2D lattice -> 'Square' or 'Hexagonal'
    :type lattice_type: string

    :param lattice_const: the distance between two particles of the lattice
    :type lattice_const: float

    :param points_file: the path to the text file with all the points where GranFilm will calculate the electric potential. |br| \
    "none" if you doesn't want the potential, "classic" if you want the potential around the particle.
    :type points_file: string

    :param materials: a list of strings for each material in media. If media = " air , mgo , ag , mgo " then materials = [ "air" , "mgo" , "ag" ]
    """

# --------------- Instance initialisation --------------- #


    def __init__(self,GRANFILM_ROOT=GRANFILM_ROOT,SOPRA_ROOT=SOPRA_ROOT,radius=8.0,truncation_ratio=0.0,radius_ratios=[1.0],
                 media="air,mgo,ag,mgo",broadening_par=0.0,broadening_perp=0.0,theta=45.0,phi=0.0,
                 energy_range=[1.5,5],pol='p',arrangement="Lattice",lattice_type="Square",
                 island_island="None",lattice_constant=20.0,multip_pos_rat=0.0,number_en_points=300,
                 multip_order=30,lower_constraint=[0.0,0.0,0.0],upper_constraint=[1e10,1.0,1e10],
                 sigma=0.005,freeze_broadening=False,points_file="None",energy=[2.2],area_ratio_pot=2.0,number_pot_points=300,
                 materials=["air","mgo","ag","mgo"],extra_parameters=[{},{},{},{}],param_file=None):
        
        self._dict = {}
        
        # Init the dictionary of the parameters

        param_dict = dict(granfilm_root=granfilm_root,sopra_root=sopra_root,radius=radius,truncation_ratio=truncation_ratio,radius_ratios=radius_ratios, media=media,
                          broadening_par=broadening_par,broadening_perp=broadening_perp,theta=theta,phi=phi, energy_range=energy_range,
                          pol=pol,arrangement=arrangement,lattice_type=lattice_type, island_island=island_island,
                          lattice_constant=lattice_constant,multip_pos_rat=multip_pos_rat,number_en_points=number_en_points, multip_order=multip_order,
                          lower_constraint=lower_constraint,upper_constraint=upper_constraint,sigma=sigma,freeze_broadening=freeze_broadening,
                          points_file=points_file,energy=energy,area_ratio_pot=area_ratio_pot,number_pot_points=number_pot_points,materials=materials,extra_parameters=extra_parameters)

        
        self.param = gf_param.GranFilm_Parameters(param_dict,param_file)
        
        # List of the data produced by the simulation

        self.list = ['info', 'epsilon', 'Susceptibilities', 'length_unit',
                     'Intensity_coeff', 'Diff_intensity_coeff', 'Polarizabilities', 'energy_unit',
                     'Fresnel', 'wavelength', 'energy','Potential']
            
        
# --------------- Dictonary-like structure --------------- #


    def __getitem__(self,index):

        assert (index in self.list), "%s is not a valid data class" % index # to make sure that pre-defined quantities can be get / set

        return self._dict[index]

    def __setitem__(self,index,value):

        assert (index in self.list), "%s is not a valid data class" % index

        self._dict[index]=value

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
        
        cmd = '%s -py -p %s -o %s' % (self.param["granfilm_root"], input_file, output_file)
        if (DEBUG=="True"): print (cmd)
        out,err,retcode,time = GranFilm_tools.run_command( cmd, tmpdir )

        # Check for errors
        if (retcode != 0):
            out,err = unicode(out, errors='replace'), unicode(err, errors='replace') # Is this needed ?
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
        if self.param["points_file"] != 'None': # Add enventual data for the potential
            data_Potential = GranFilm_tools.get_results_from_file( os.path.join(tmpdir,"Results_potential.dat")  )
        
        # Clean up the un-needed directory structure
        shutil.rmtree(tmpdir)   
                
        # Populate the data container   
        GranFilm_results.GranFilm_Results(self, data, data_RefTrans, data_Epsilon, data_Potential)
        
        # Print some final comments
        if (DEBUG=="True"): print ("GranFilm dictionary populated with simulation results!")

        return
    

# --------------- Other methods --------------- #

    
    def plot(self,coeff_type,param,abscisse="energy",legend=""):

        x=self[abscisse]
        y=self[coeff_type](*param)

        if legend=="":

            plt.plot(x,y)

        else:

            plt.plot(x,y,label=legend)
            plt.legend()

        return
