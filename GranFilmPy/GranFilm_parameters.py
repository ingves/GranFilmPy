import numpy as np
import os
import ast
import matplotlib.pyplot as plt

python_interface = os.path.dirname(os.path.realpath(__file__)).replace('\\','/') # Path to the directory of the .py files constituting the python interface for GranFilm   
    
class GranFilm_Parameters(dict):


# --------------- Init method --------------- #


    def __init__(self,param_dict,potential_dict,param_file):

        # --- Initialisation of the gf.param dictionary

        if param_file is None: # Take the parameters from the param_dict
            
            self.update(param_dict)
            self.Potential = GranFilm_Parameters_Potential(potential_dict)
        
        else: # Take parameters from the param_file and complain if the file is not correct

            self.update(param_dict)
            self.Potential = GranFilm_Parameters_Potential(potential_dict)
            self.File(param_file,'r')            
                
    def __getitem__(self,key):

        assert (key in self.keys()), "%s is not a valid parameter" % key # To make sure that only pre-defined keys can be get / set
        return self.get(key)

    def __setitem__(self,key,value):

        assert (key in self.keys()), "%s is not a valid parameter" % key
        dict.__setitem__(self,key,value)
        
    def __str__(self,format="Default"):
        
       
        if format == "Default":
            
            ret = "\n----- Standard parameters -----\n\n"
            
            ret += "&Global \n"
            ret += "  GRANFILM_ROOT     = '{}' \n".format(self["GRANFILM_ROOT"])
            ret += "  SOPRA_ROOT        = '{}' \n".format(self["SOPRA_ROOT"])
            ret += "/\n\n"
            
            ret += "&Source \n"
            ret += "  theta             = {} \n".format(self["theta"])
            ret += "  phi               = {} \n".format(self["phi"])
            ret += "  pol               = '{}' \n".format(self["pol"])
            ret += "  energy_range      = {} \n".format(self["energy_range"])
            ret += "/\n\n"

            ret += "&Geometry \n"
            ret += '  radius            = {} \n'.format(self["radius"])
            ret += '  truncation_ratio  = {} \n'.format(self["truncation_ratio"])
            ret += '  radius_ratios     = {} \n'.format(self["radius_ratios"])
            ret += '  media             = {} \n'.format(self["media"])
            ret += "/\n\n"

            ret += "&Interaction \n"
            ret += "  arrangement       = '{}' \n".format(self["arrangement"])
            ret += "  lattice_type      = '{}' \n".format(self["lattice_type"])
            ret += "  island_island     = '{}' \n".format(self["island_island"])
            ret += "  lattice_constant  = {} \n".format(self["lattice_constant"])
            ret += "/\n\n"
            
            ret += "&Numerics \n"
            ret += "  multip_pos_rat    = {} \n".format(self["multip_pos_rat"])
            ret += "  number_en_points  = {} \n".format(self["number_en_points"])
            ret += "  multip_order      = {} \n".format(self["multip_order"])
            ret += "/\n\n"
            
            ret += "----- Parameters for potential calculation -----\n\n"

            ret += "&Potential \n"
            ret += "  points_file       = '{}' \n".format(self.Potential["points_file"])
            ret += "  energy            = {} \n".format(self.Potential["energy"])
            ret += "  area_ratio_pot    = {} \n".format(self.Potential["area_ratio_pot"])
            ret += "  number_pot_points = {} \n".format(self.Potential["number_pot_points"])
            ret += "/"            

        elif format == "Fortran":
            
            ret  = "&Global \n" # Global
            ret += " SOPRA_ROOT = '%s' \n" % self["SOPRA_ROOT"]
            ret += "/\n\n"

            ret += "&Source \n"
            ret += "  Theta0           = %f \n"    % self["theta"]
            ret += "  Phi0             = %f \n"    % self["phi"]
            ret += "  Polarization     = '%s' \n"  % self["pol"]
            ret += "  Energy_Range     = "+str(self["energy_range"][0])+", "+str(self["energy_range"][1])+"\n"
            ret += "/\n\n"

            ret += "&Geometry \n" # Geometry          
            # Change the radius into a string and remove potential brackets so that 8.0 -> "8.0" and [8.0,5.0] -> "8.0,5.0"
            try:
                ret += '  Radius           = "{},{}" \n'.format(self["radius"][1],self["radius"][0])
            except:
                ret += '  Radius           = "{}" \n'.format(self["radius"])                
            ret += "  Truncation_Ratio = %f \n"      % self["truncation_ratio"]
            
            # Change the radius_ratios into a string and remove potential brackets
            ret += '  Radius_Ratios    = "{}" \n'.format(str(self["radius_ratios"]).replace('[','').replace(']',''))                                                            
            ret += '  Media            = "%s" \n'    % ','.join(self["media"]) # Convert the list of media into a string
            ret += "/\n\n"

            ret += "&Interaction \n" # Interaction
            ret += "  Arrangement               = '%s' \n"    % self["arrangement"]
            ret += "  Lattice_Type              = '%s' \n"    % self["lattice_type"]
            ret += "  Island_Island_Interaction = '%s' \n"    % self["island_island"]
            ret += "  Lattice_Constant          = %f \n"      % self["lattice_constant"]
            ret += "/\n\n"

            ret += "&Numerics \n" # Numerics
            ret += "  Multipole_Position_Ratio  = %f \n"    % self["multip_pos_rat"]
            ret += "  No_Energy_Points          = %i \n"    % self["number_en_points"]
            ret += "  Multipole_Order           = %i \n"    % self["multip_order"]
            ret += "/\n\n"

            if self.Potential["points_file"] != "None": # Potential
                ret += "&Potential \n"
                if self.Potential["points_file"] in ["classic","surface"]:
                    ret += "  Points_File = '{}' \n".format(python_interface+'/'+self.Potential["points_file"])
                else:
                    ret += "  Points_File = '%s' \n"   % os.path.abspath(self.Potential["points_file"])
                ret += "  Energy = %s \n" % ', '.join([str(x) for x in self.Potential["energy"]])
                ret += "/\n\n"
                
#            for mat in self["media"].split(','):
#                ret += "&%s \n"     % mat
#                try :
#                    ret += "  Material = \"%s\" \n"    % mat.lower()
#                except:
#                    ret += "  Material = \"%s\" \n"    % mat.toLower()
#                finally:
#                    for key in self["extra_parameters"][k].keys():
#                        if key == "Epsilon_Scale":
#                            ret += "  {} = {} \n".format(key,str(self["extra_parameters"][k][key]).replace('[','').replace(']',''))
#                        else:
#                            ret += "  {} = {} \n".format(key,self["extra_parameters"][k][key])
#            ret += "/\n\n"
        else:
            ret = "Wrong format for parameters printing. Try format = 'Default' or 'Fortran'."
        return ret

    def __call__(self,format="Default"):

        print(self.__str__(format))


# --------------- File method to read / write an external input file --------------- #


    def File(self, fname, mode="None",format="Extended"):
        """
        Reads/Writes parameters to/from an external file that is in the sif-format. 

        Options
        -------
        fname : string
           Name of the file on disk. The file is in the sif-format
        mode : string 
           Mode specifies file reading ('r') or writing ('w') 
        
        """

        # --- Writing parameter file
        
        if (mode=="w"):
            try :
                with open(fname,'w') as ofile:
                    ofile.write(str(self.__str__(format)))
            except:
                print (" ERROR : Writing file %s "% fname)
            return    

        # --- Reading parameter file
        
        elif (mode=="r"):

            try:
                
                assert(os.path.isfile(fname)), "ERROR: param_file = {} does not exist".format(fname)

                with open(fname,'r') as f:
                    L = f.readlines()
                    
                L2 = []

                for s in L:
                    if not(s[0] in ['&','/'] or s=='\n'): # remove &Geometry, \n and /
                        L2.append(s)

                L3=[s.replace('\n','').replace(' ','').split('=') for s in L2] # remove blanks and split in (param,value)
                D = dict(L3)
                param_keys = self.keys()
                potential_keys = self.Potential.keys()
                
                for key in D.keys():
                    
                    if key in param_keys: # To only update existing parameters
                        self[key] = ast.literal_eval(D[key]) # turn "[1,2,3]" to [1,2,3] and is safe
                    elif key in potential_keys:
                        self.Potential[key] = ast.literal_eval(D[key])
                    else:
                        print("Warning: {} is not a valid parameter and will not be exported.\n".format(key))
                
            except ValueError:

                print("""ERROR: The format of {} is not correct for the export of the parameters.\n
                      The format is: parameter1=value1\nparameter2=value2\n...\n\n
                      Parameters are unchanged.\n""".format(fname))

        else:
            print("ERROR: mode = {} must be 'w' or 'r'".format(mode))


    def Material(self,material,x_axis="eV"):
        
        path = self["SOPRA_ROOT"] + material + ".nk" # path to the .nk file of the medium
        
        try:
            assert(x_axis in ['eV','microns']), "ERROR: x_axis = {} is not valid. Try x_axis = 'eV' or 'microns'".format(x_axis)
            with open(path,'r') as f:
                s = [ast.literal_eval(v) for v in f.readline().split()] # gets the format of the file (1 for eV, 2 for microns), the interval and the number of points
                x = np.linspace(s[1],s[2],s[3]+1)               
                y = []        
                for k in range(s[3]+1):     # fills in list with values of the refractive index
                    a = f.readline().split()
                    y.append(float(a[0])+1j*float(a[1]))
            
            ref_index = np.array(y,dtype=complex)
            epsilon = ref_index**2
            
            if not((x_axis == 'microns' and s[0] == 2) or (x_axis == 'eV' and s[0] == 1)):
                x = 1.24/x # if wrong unit for x, then convert x to eV or microns
            
            plt.subplot(121) # Plot of the dielectric function
            plt.plot(x,epsilon.real,label="Real part")
            plt.title("Dielectric function of %s" %material)
            plt.plot(x,epsilon.imag,label="Imag part")
            if x_axis == "eV":
                plt.xlabel("Energy (eV)")
            else:
                plt.xlabel("Wavelength (microns)")
            plt.legend()
            
            plt.subplot(122) # Plot of the refractive index
            plt.plot(x,ref_index.real,label="Real part")
            plt.title("Refractive index of %s" %material)
            plt.plot(x,ref_index.imag,label="Imag part")
            if x_axis == "eV":
                plt.xlabel("Energy (eV)")
            else:
                plt.xlabel("Wavelength (microns)")
            plt.legend()
            
            plt.show()

        except:
            print("ERROR: Reading file {}".format(path))

    def write_points_file(self):

        # --- Initialisation of the input file for potential

        n_pts = self.Potential["number_pot_points"] # number of points for the potential
        area_ratio = self.Potential["area_ratio_pot"]
        radius = self["radius"]
        points_file = self.Potential["points_file"]
        
        try:
            x_radius = radius[0]
            y_radius = radius[1]

        except:
            x_radius = radius
            y_radius = radius

        if points_file == "classic": # Potential calculated around the particle

            potential_file = open(python_interface+"/classic","w")
            string = ""
            
            # Set the spatial area (rectangle around the particle)
            
            x = np.linspace(-area_ratio*x_radius,area_ratio*x_radius,n_pts)
            y = np.linspace(-area_ratio*y_radius,area_ratio*y_radius,n_pts)

            for v1 in x:
                for v2 in y:
                    string += "{} 0. {}\n".format(v1,v2)
                    
            potential_file.write(string)      
            potential_file.close()

        elif points_file == "surface": # Potential calculated at the surfaces of the coatings

            string = ""
            potential_file = open(python_interface+"/surface","w")
            theta = np.linspace(0,2*np.pi,n_pts)

            # Set the spatial area (here all the coatings surfaces)
            
            for c in self["radius_ratios"]:
                for t in theta:
                    string += "{} 0. {}\n".format((0.001 + x_radius)*np.cos(t),-(0.001 + y_radius)*np.sin(t))
                    string += "{} 0. {}\n".format((x_radius-0.001)*np.cos(t),-(y_radius-0.001)*np.sin(t))

            potential_file.write(string)
            potential_file.close()

        elif (points_file != 'None') and (not os.path.isfile(points_file)): # points_file doesn't exist

            print("points_file = %s does not exist.\nPotential will not be calculated (points_file = 'None')" % self.Potential["points_file"])
            self.Potential["points_file"] = 'None'

class GranFilm_Parameters_Potential(dict):
    
    """ 
        Class which handle the parameters for the calculation of the potential.
        Works like the GranFilm_Parameters class.
        
    """
    
    def __init__(self,potential_dict):
        
        self.update(potential_dict)
        
    def __getitem__(self,key):

        assert (key in self.keys()), "%s is not a valid parameter" % key # To make sure that only pre-defined keys can be get / set
        return self.get(key)

    def __setitem__(self,key,value):

        assert (key in self.keys()), "%s is not a valid parameter" % key
        dict.__setitem__(self,key,value)
