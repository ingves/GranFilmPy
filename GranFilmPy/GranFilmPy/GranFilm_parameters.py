import numpy as np
import os
import ast
import matplotlib.pyplot as plt

python_interface = os.path.dirname(os.path.realpath(__file__)).replace('\\','/') # Path to the directory of the .py files constituting the python interface for GranFilm   
    
class GranFilm_Parameters(dict):


# --------------- Special methods --------------- #


    def __init__(self,param_dict,param_file):

        # --- Initialisation of the gf.param dictionary

        if param_file is None: # Take the parameters from the param_dict
            
            self.update(param_dict)
            
        else: # Take parameters from the param_file and complain if the file is not correct

            self.update(param_dict)
            self.File(param_file,'r')
            
                
    def __getitem__(self,key):

        assert (key in self.keys()), "%s is not a valid parameter" % key # To make sure that only pre-defined keys can be get / set

        return self.get(key)

    def __setitem__(self,key,value):

        assert (key in self.keys()), "%s is not a valid parameter" % key

        dict.__setitem__(self,key,value)
        
    def __str__(self,format="Python"):
        
        if format == "Python":
            
            ret = "&Global \n"
            ret += "  granfilm_root     = '{}' \n".format(self["granfilm_root"])
            ret += "  sopra_root        = '{}' \n".format(self["sopra_root"])
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
            ret += '  broadening_par    = {} \n'.format(self["broadening_par"])
            ret += '  broadening_perp   = {} \n'.format(self["broadening_perp"])
            ret += '  radius_ratios     = {} \n'.format(self["radius_ratios"])
            ret += '  media             = "{}" \n'.format(self["media"])
            ret += "/\n\n"

            ret += "&Interaction \n"
            ret += "  arrangement       = '{}' \n".format(self["arrangement"])
            ret += "  lattice_type      = '{}' \n".format(self["lattice_type"])
            ret += "  island_island     = '{}' \n".format(self["island_island"])
            ret += "  lattice_constant  = {} \n".format(self["lattice_constant"])
            ret += "/\n\n"

            ret += "&Curvefitting \n"
            ret += "  lower_constraint  = {} \n".format(self["lower_constraint"])
            ret += "  upper_constraint  = {} \n".format(self["upper_constraint"])
            ret += "  sigma             = {} \n".format(self["sigma"])
            ret += "  freeze_broadening = {} \n".format(self["freeze_broadening"])
            ret += "/\n\n"

            ret += "&Potential \n"
            ret += "  points_file       = '{}' \n".format(self["points_file"])
            ret += "  energy            = {} \n".format(self["energy"])
            ret += "  area_ratio_pot    = {} \n".format(self["area_ratio_pot"])
            ret += "  number_pot_points = {} \n".format(self["number_pot_points"])
            ret += "/\n\n"

            ret += "&Numerics \n"
            ret += "  multip_pos_rat    = {} \n".format(self["multip_pos_rat"])
            ret += "  number_en_points  = {} \n".format(self["number_en_points"])
            ret += "  multip_order      = {} \n".format(self["multip_order"])
            ret += "/\n\n"

            ret += "&Media \n"
            ret += "  materials         = {} \n".format(self["materials"])
            ret += "  extra_parameters  = {} \n".format(self["extra_parameters"])
            ret += "/"

        else:
            
            ret  = "&Global \n" # Global
            ret += " SOPRA_ROOT = '%s' \n" % self["sopra_root"]
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
            ret += "  Broadening_Par   = %f \n"      % self["broadening_par"]
            ret += "  Broadening_Perp  = %f \n"      % self["broadening_perp"]
            # Change the radius_ratios into a string and remove potential brackets
            ret += '  Radius_Ratios    = "{}" \n'.format(str(self["radius_ratios"]).replace('[','').replace(']',''))                                                            
            ret += '  Media            = "%s" \n'    % self["media"]
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

            ret += "&Curvefitting \n" # Curvefitting
            ret += ("  Lower_Constraint = %s \n"   % 
                            ', '.join([str(x) for x in self["lower_constraint"]]))
            ret += ("  Upper_Constraint  = %s \n"  % 
                            ', '.join([str(x) for x in self["upper_constraint"]]))
            ret += "  sigma  = %f \n" % self["sigma"]
            ret += "  freeze_broadening  = {} \n".format(self["freeze_broadening"])
            ret += "/\n\n"

            if self["points_file"] != "None": # Potential
                ret += "&Potential \n"
                if self["points_file"] in ["classic","surface"]:
                    ret += "  Points_File = '{}' \n".format(python_interface+'/'+self["points_file"])
                else:
                    ret += "  Points_File = '%s' \n"   % os.path.abspath(self["points_file"])
                ret += "  Energy = %s \n" % ', '.join([str(x) for x in self["energy"]])
                ret += "/\n\n"
                
            for k in range(len(self["materials"])): # Materials
                mat = self["materials"][k]
                ret += "&%s \n"     % mat
                try :
                    ret += "  Material = \"%s\" \n"    % mat.lower()
                except:
                    ret += "  Material = \"%s\" \n"    % mat.toLower()
                finally:
                    for key in self["extra_parameters"][k].keys():
                        if key == "Epsilon_Scale":
                            ret += "  {} = {} \n".format(key,str(self["extra_parameters"][k][key]).replace('[','').replace(']',''))
                        else:
                            ret += "  {} = {} \n".format(key,self["extra_parameters"][k][key])
            ret += "/\n\n"

        return ret

    def __call__(self,format="Python"):

        print(self.__str__(format))


# --------------- Other methods --------------- #


    def File(self, fname, mode="None",format="Python"):
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
                D,param_dict = dict(L3),{}
                keys = self.keys()
                
                for key in D.keys():
                    if key in keys: # To only update existing parameters
                        param_dict[key] = ast.literal_eval(D[key]) # turn "[1,2,3]" to [1,2,3] and is safe
                    else:
                        print("Warning: {} is not a valid parameter and will not be exported.\n".format(key))

                self.update(param_dict)
                
            except ValueError:

                print("ERROR: The format of {} is not correct for the export of the parameters.\nGood format is: parameter=value\nParameters are unchanged.\n".format(fname))

        else:
            print("ERROR: mode = {} must be 'w' or 'r'".format(mode))


    def Epsilon(self,material):
        
        path = self["sopra_root"] + material + ".nk"
        
        try:
            with open(path,'r') as f:

                s =[float(v) for v in f.readline().split()[1:]]
                print("Energy_min = {}, Energy_max = {}".format(s[0],s[1]))
                x = np.linspace(s[0],s[1],s[2]+1)
                L = f.readlines()
                
            ref_index = []
        
            for line in L:
                a = line.split()
                ref_index.append(float(a[0])+1j*float(a[1]))

            ref_index = np.array(ref_index,dtype=complex)
            epsilon = ref_index**2
            
            plt.subplot(211)
            plt.plot(x,epsilon.real,label="Real")
            plt.title("Dielectric function of %s" %material)
            plt.legend()
            plt.subplot(212)
            plt.plot(x,epsilon.imag,label="Imag",color='r')
            plt.xlabel("Energy (eV)")
            plt.legend()
            plt.show()

        except:
            print("ERROR: Reading file {}".format(path))

    def write_points_file(self):

        # --- Initialisation of the input file for potential

        n_pts = self["number_pot_points"] # number of points for the potential
        area_ratio = self["area_ratio_pot"]
        radius = self["radius"]
        points_file = self["points_file"]
        
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

            print("points_file = %s does not exist.\nPotential will not be calculated (points_file = 'None')" % self["points_file"])
            self["points_file"] = 'None'
        
##    def check_parameters(self):
##
##        assert(os.path.isfile(self["granfilm_root"])),"ERROR: The directory %s does not exist.\n" %self["granfilm_root"]
##       
##        assert(os.path.isdir(self["sopra_root"])),"ERROR: The file %s does not exist.\n" %self["sopra_root"]
##       
##        assert(0<=self["theta"]<90),"ERROR: theta = {} must satisfy 0 <= theta < 90 degrees.\n" %self["theta"]
##       
##        assert(self["pol"] in ['p','s']),"ERROR: polarization = {} must be 'p' or 's'.\n".format(self["pol"])
##       
##        L1,L2 = [],[]
##        for mat in self["materials"]:
##            path = self["sopra_root"]+mat+".nk"
##            assert(os.path.isfile(path)),"ERROR: %s does not exist.\n" %path
##            f = open(path,'r')
##            s = f.readlines()[0].split()
##            L1.append(s[1]) # min
##            L2.append[s[2]) # max
##
##        energy_min,energy_max = max(L1),min(L2)        
##        en_range = self["energy_range"]
##        assert(energy_min<=en_range[0] and en_range[1]<=energy_max),"ERROR: energy_range = {} must be included in [{},{}].\n".format(en_range,energy_min,energy_max)
##
##        r = self["radius"]
##        assert(0<r<=self["lattice_constant"]),"ERROR: radius = {} must satisfy 0 < radius <= lattice_constant = {}.\n".format(self["rad
##
        
