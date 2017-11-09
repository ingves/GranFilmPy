import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc,Ellipse
from GranFilm_tools import reformat
from matplotlib.colors import LinearSegmentedColormap
from skimage import measure


# --------------- Function to populate the data structure --------------- #


def GranFilm_Results(gf,data,data_RefTrans,data_Epsilon,data_Potential):

    import GranFilm_tools

    # Set the info attribute
    gf["info"]        =  " Time of generation : %s " % GranFilm_tools.timestamp() 

    # Units
    gf["length_unit"] = "nano meters"   # Update this 
    gf["energy_unit"] = "eV"            # Update this 
    
    # Set energy and wavelength
    gf["energy"]     = data[:,0]
    gf["wavelength"] = data[:,19]

    # Fresnel coefficients
    gf["Fresnel"] = GranFilm_Results_Fresnel( data_RefTrans )

    # Epsilon
    gf["epsilon"] = GranFilm_Results_Epsilon( data_Epsilon, gf )

    # Suceptibilities
    gf["Susceptibilities"] =  GranFilm_Results_Suceptibilities( data )

    # Polarizabilities
    gf["Polarizabilities"] =  GranFilm_Results_Polarizabilities( data )

    # Intensity coefficients

    gf["Intensity_coeff"] = GranFilm_Results_Intensity_Coeff(gf)

    # Differential Intensity coefficients
    media = gf.param["media"].split(',')
    if media[0] != media[1]:
        gf["Diff_intensity_coeff"]  = GranFilm_Results_Diff_Intensity_Coeff(gf)
    
    # Potential
    if data_Potential is not None:
        gf["Potential"]  = GranFilm_Results_Potential(gf,data_Potential)

    
# --------------- Diff_Intensity_Coeff --------------- #

    
class GranFilm_Results_Diff_Intensity_Coeff:

    def __init__(self,gf):

        I_coeff = gf["Intensity_coeff"]
        
        self.dRp_Rp   = (I_coeff('R','p') - I_coeff('R','p','Flat'))/I_coeff('R','p','Flat')
        self.dRs_Rs   = (I_coeff('R','s') - I_coeff('R','s','Flat'))/I_coeff('R','s','Flat')
        self.dRps_Rps = (I_coeff('R','ps') - I_coeff('R','ps','Flat'))/I_coeff('R','ps','Flat')

        self.dTp_Tp   = (I_coeff('T','p') - I_coeff('T','p','Flat'))/I_coeff('T','p','Flat')
        self.dTs_Ts   = (I_coeff('T','s') - I_coeff('T','s','Flat'))/I_coeff('T','s','Flat')
        self.dTps_Tps = (I_coeff('T','ps') - I_coeff('T','ps','Flat'))/I_coeff('T','ps','Flat')

        self.dAp_Ap   = (I_coeff('A','p') - I_coeff('A','p','Flat'))/I_coeff('A','p','Flat')
        self.dAs_As   = (I_coeff('A','s') - I_coeff('A','s','Flat'))/I_coeff('A','s','Flat')
        self.dAps_Aps = (I_coeff('A','ps') - I_coeff('A','ps','Flat'))/I_coeff('A','ps','Flat')


    def __call__(self,ref_trans,pol):
        
        # --- some assertions
        assert (ref_trans in {'R','T'}), "ERROR : Unsupported coefficient -> Valid options: 'R' , 'T'"
        assert (pol in {'p','s','ps'}), "ERROR : Unsupported polarisation -> Valid options: 'p' , 's' , 'ps'"
        
        return getattr(self,'d'+ref_trans+pol+'_'+ref_trans+pol)


# --------------- Intensity_Coeff --------------- #

    
class GranFilm_Results_Intensity_Coeff:
    

    def __init__(self,gf):

        fresnel = gf["Fresnel"]
        theta = gf.param["theta"]*np.pi/180.0
        media=gf.param["media"].split(',')
        eps_air = gf["epsilon"](media[0])
        eps_substrate = gf["epsilon"](media[1])

        # ----- Reflection ----- #
        rp=fresnel('r','p')
        rp_flat=fresnel('r','p','Flat')
        rs=fresnel('r','s')
        rs_flat=fresnel('r','s','Flat')

        self.Rp = np.abs(rp)**2.0
        self.Rp_flat = np.abs(rp_flat)**2.0
        self.Rs = np.abs(rs)**2.0
        self.Rs_flat = np.abs(rs_flat)**2.0
        self.Rps = (self.Rp + self.Rs)/2.0
        self.Rps_flat = (self.Rp_flat + self.Rs_flat)/2.0
        
        # ----- Transmission ----- #
        tp=fresnel('t','p')
        tp_flat=fresnel('t','p','Flat')
        ts=fresnel('t','s')
        ts_flat=fresnel('t','s','Flat')
        e = np.abs((eps_substrate/eps_air)**0.5)
        a = np.cos(np.arcsin(np.sin(theta)/e))/np.cos(theta)
        self.Tp = e*a*np.abs(tp)**2.0
        self.Tp_flat = np.abs(tp_flat)**2.0
        self.Ts = e*a*np.abs(ts)**2.0
        self.Ts_flat = np.abs(ts_flat)**2.0
        self.Tps = (self.Tp + self.Ts)/2.0
        self.Tps_flat = (self.Tp_flat + self.Ts_flat)/2.0

        # ----- Absorption ----- #
        self.Ap = 1-self.Rp-self.Tp
        self.Ap_flat = 1-self.Rp_flat-self.Tp_flat
        self.As = 1-self.Rs-self.Ts
        self.As_flat = 1-self.Rs_flat-self.Ts_flat
        self.Aps = 1-self.Rps-self.Tps
        self.Aps_flat = 1-self.Rps_flat-self.Tps_flat
        

    def __call__(self,ref_trans,pol,geometry='Island'):

        # --- some assertions
        assert (ref_trans in {'A','R','T'}), "ERROR : Unsupported coefficient -> Valid options: 'R' , 'T'"
        assert (pol in {'p','s','ps'}), "ERROR : Unsupported polarisation -> Valid options: 'p' , 's' , 'ps'"  
        assert (geometry in {'Island','Flat'}), "ERROR : Unsupported geometry -> Valid options: 'Island' , 'Flat'"

        if geometry=='Island':
            return getattr(self,ref_trans+pol)
        else:
            return getattr(self,ref_trans+pol+'_flat')


# --------------- Fresnel_Coeff --------------- #


class GranFilm_Results_Fresnel:
    
    """
    Class for populating the data container

    Attributes:

    Methods:
        None at this point.
  
    Use:
        Call the class to read data

    """
    
    def __init__(self, data):
        
        imu = 1j                         # imaginar unit
        #
        self.energy = data[:,0]
        # reflction amplitudes
        self.r_p       = data[:,1]  + imu * data[:,2]
        self.r_p_flat  = data[:,3]  + imu * data[:,4]
        self.r_s       = data[:,5]  + imu * data[:,6]
        self.r_s_flat  = data[:,7]  + imu * data[:,8]
        # transmission amplitudes
        self.t_p       = data[:,9]  + imu * data[:,10]
        self.t_p_flat  = data[:,11] + imu * data[:,12]
        self.t_s       = data[:,13] + imu * data[:,14]
        self.t_s_flat  = data[:,15] + imu * data[:,16]
        

    def __call__(self,ref_trans,polarization,geometry="Island",fmt="Complex",deg="True"): 

        """
        Gets the fresnel amplitude for reflection or transmission for a given linear polarization 
        of the incident light           

        Options
        -------
        ref_trans : string
            Reflection ('r') or transmission ('t') amplitudes
        polarization : string
            P-polarization ('p') or s-polarization ('s') of the incident light.
            If polarization='ps' the sum of p and s-polarization is returned
        geometry : string
            By default ('Island'), the amplitude for the supported granular system is returned.
            If geometry="Flat", the classic Fresnel amplitudes for the corresponding flat geometry is returned. 
        
        """
        
        # --- some assertions
        assert (ref_trans in {'r','t'}), "ERROR : Unsupported coefficient -> Valid options: 'R' , 'T'"
        assert (polarization in {'p','s','ps'}), "ERROR : Unsupported polarization -> Valid options: 'p' , 's' , 'ps'"  
        assert (geometry in {'Island','Flat'}), "ERROR : Unsupported geometry -> Valid options: 'Island' , 'Flat'"


        # ----- Reflection ----- #
        
        if (ref_trans=="r"):
            # Initial value
            Fresnel_Amplitude = 0.   
            # --- p-polarization
            if (polarization in {'ps','p'}):
                if (geometry=="Flat"):
                    Fresnel_Amplitude =+  self.r_p_flat
                else:
                    Fresnel_Amplitude =+  self.r_p
            # --- s-polarization
            if (polarization in {'ps','s'}):
                if (geometry=="Flat"):
                    Fresnel_Amplitude =+  self.r_s_flat
                else:
                    Fresnel_Amplitude =+  self.r_s
            # Return result
            return reformat(Fresnel_Amplitude,fmt=fmt,deg=deg)
        
        # ----- Transmission ----- #
        
        if (ref_trans=="t"):
            # Initial value
            Fresnel_Amplitude = 0.   
            # --- p-polarization
            if (polarization in {'ps','p'}):
                if (geometry=="Flat"):
                    Fresnel_Amplitude =+  self.t_p_flat
                else:
                    Fresnel_Amplitude =+  self.t_p
            # --- s-polarization
            if (polarization in {'ps','s'}):
                if (geometry=="Flat"):
                    Fresnel_Amplitude =+  self.t_s_flat
                else:
                    Fresnel_Amplitude =+  self.t_s
            # Return result
            return Fresnel_Amplitude


# --------------- Susceptibilities --------------- #


class GranFilm_Results_Suceptibilities:
    
    """
    Class for populating the data container

    Attributes:

    Methods:
        Not used at this point.

    Use:
         Call the class to read data

    """
    
    def __init__(self, data):
                
        imu = 1j                         # imaginar unit
        #
        self.energy = data[:,0]
        #
        self. gamma   = data[:,11]  + imu * data[:,12]
        self. beta    = data[:,13]  + imu * data[:,14]
        #
        self. delta   = data[:,15]  + imu * data[:,16]
        self. tau     = data[:,17]  + imu * data[:,18]


    def __call__(self, type="None", fmt='Complex', deg='True'):
        """
        Get the susceptibility of the system

        Options
        -------
        type: string
            which surface susceptibility to return (delta, bets,gamma or tau)
            Default : None
        fmt : string
           Format of the returned value {'Complex', 'Real', 'Imag', 'Amplitude', 'Phase' }
           Default : Complex
        deg: string
           If fmt=Phase use deg (or not) 
        """           

        # error checking
        FMT={'Complex', 'Real', 'Imag', 'Amplitude', 'Phase' }
        assert (type in {'delta','beta', 'gamma', 'tau'}), "ERROR : type=%s is not supported (valid options: 'Dipole','Quadrupole') " % type        
        assert (fmt in FMT), "ERROR : type=%s is not supported " % fmt

        if (type=='delta'):
            susceptibility = self.delta
        elif (type=='beta'):
            susceptibility = self.beta
        elif (type=='gamma'):
            susceptibility = self.gamma
        elif (type=='tau'):
            susceptibility = self.tau

        # Output format
        return reformat(susceptibility, fmt=fmt, deg=deg)


# --------------- Polarizabilities --------------- #


class GranFilm_Results_Polarizabilities:
    
    """
    Class for populating the data container

    Attributes:

    Methods:
       Not used at this point.

    Use:
        Call the class to read data

    """
    
    def __init__(self, data):
                
        imu = 1j                         # imaginar unit
        #
        self.energy = data[:,0]
        #
        self.alpha_parallel       = data[:,3]  + imu * data[:,4]
        self.alpha_perpendicular  = data[:,5]  + imu * data[:,6]
        #
        self.alpha_10_parallel       = data[:,9]   + imu * data[:,10]
        self.alpha_10_perpendicular  = data[:,11]  + imu * data[:,12]


    def __call__(self, orientation="None", type='Dipole', fmt='Complex', deg='True'):
        """
        Get the polarizability of the system

        Options
        -------
        orientation : string
           'Parallel' or 'Perpendicular' to the surface of the substrate
        type: string
            Polarizability type to return ('Dipole' or 'Quadrupole')
            Default : Dipole
        fmt : string
           Format of the returned value {'Complex', 'Real', 'Imag', 'Amplitude', 'Phase' }
           Default : Complex
        deg: string
           If fmt=Phase use deg (or not) 
        
        """                

        # make sure that the input have supported values
        FMT={'Complex', 'Real', 'Imag', 'Amplitude', 'Phase' }
        assert (orientation in {'Parallel','Perpendicular'}), "ERROR : type=%s is not supported (valid options: 'Parallel','Perpendicular') " % orientation 
        assert (type in {'Dipole','Quadrupole'}), "ERROR : type=%s is not supported (valid options: 'Dipole','Quadrupole') " % type        
        assert (fmt in FMT), "ERROR : type=%s is not supported (valid options: 'Complex','Real','Imag','Amplitude','Phase'" % fmt


        # --- Get data
        # 
        # Parallel Polarizability
        if (orientation=='Parallel'):
            if (type=='Dipole'):     polarizability = self.alpha_parallel 
            if (type=='Quadrupole'): polarizability = self.alpha_10_parallel 
        #
        # Parallel Polarizability
        if (orientation=='Perpendicular'):
            if (type=='Dipole'):     polarizability = self.alpha_perpendicular
            if (type=='Quadrupole'): polarizability = self.alpha_10_perpendicular

            
        # Output format
        return reformat(polarizability, fmt=fmt, deg=deg)



# --------------- Dielectric functions of the materials --------------- #


class GranFilm_Results_Epsilon(dict):
    
    """
    Class for populating the data container

    Attributes:

    Methods:
        Not used at this point.

    Use:
        Call the class to read data

    """
    
    def __init__(self, data, gf):

        # Imaginary unit
        imu = 1j 

        # --- loop over materials
        i=-1
        for material in gf.param["media"].split(','):
            i += 2
            self[material] = data[:,i]  + imu * data[:,i+1]


    def __call__(self,medium="None"):
        """
        Get the dielectric function of 'medium' that must be in the list defined media 
        
        Options
        -------
        medium : string
            Medium id for one of the media in the geometry
   
        """        
        if ( medium in self.keys() ):
            return self[medium]
        else:
            print ("ERROR : Defined materials are : %s " % self.keys())
            

# --------------- Potential --------------- #

class GranFilm_Results_Potential():

    def __init__(self,gf,data_Potential):

        self.data = data_Potential
        self.points_file = gf.param["points_file"]
        
        if self.points_file == "classic":
            
            x_border = self.data[0][0] # -area_ratio * x_radius
            y_border = self.data[0][2] # -area_ratio * y_radius
            self.x_radius = x_border/-gf.param["area_ratio_pot"]
            self.y_radius = y_border/-gf.param["area_ratio_pot"]
            self.area = (x_border,-x_border,-y_border,y_border)
            self.truncation_ratio = gf.param["truncation_ratio"]
            self.interface_level = self.truncation_ratio*self.y_radius
            self.media = gf.param["media"].split(',')
            self.radius_ratios = gf.param["radius_ratios"]

            # Reshapment of the data to create the "map" with the values of the potential
            
            n=len(self.data)
            p=int(np.sqrt(n))
            P=np.zeros((p,p),dtype=complex)
            
            for i in range(p):
                for k in range(p):
                    P[k,i]=self.data[i*p+k][4]+1j*self.data[i*p+k][5]
                    
            self.potential = P

        elif self.points_file == "surface":
            
            # Calculation of the error
            s = 0
            n = len(self.data)
            for k in range(n//2):
                s += np.abs(self.data[2*k][4]-self.data[2*k+1][4])

            self.error = 2*s/n
            
        
    def __call__(self,output='plot',fmt="Real",deg='True',equipotential=False):

        
        if output=='plot': # Plot the potential

            formated_potential = reformat(self.potential,fmt=fmt,deg=deg)

            assert (self.points_file == "classic"), "ERROR: Potential can only be plotted for points_file = classic"
            fig = plt.figure()
            ax = fig.add_subplot(111)
            x_max = self.area[1]
            y_max = self.area[3]
            
            plt.plot((-x_max,x_max),(self.interface_level,self.interface_level),color='k',linestyle='dashed',linewidth=1.0) # horizontal line
            ax.text(0,0,'.',fontsize=20,clip_on=True) # center of the particle
            plt.xlabel("x (nm)")
            plt.ylabel("y (nm)")
            plt.title("Electric potential around the particle")

            im = plt.imshow(formated_potential,cmap=plt.cm.RdBu,extent=self.area) # plot the potential map
            plt.colorbar(im)

            # add coatings and text for media
            ax.text(-0.9*x_max,0.9*y_max,self.media[0],clip_on=True) # medium 1
            ax.text(-0.9*x_max,-0.9*y_max,self.media[1],clip_on=True) # medium 2
            A=[]
            B=[]
            
            for k in range(len(self.radius_ratios)): # plot the different coatings
                c=self.radius_ratios[k]
                plot_coating(self.x_radius*c,self.y_radius*c,self.media[2*k:2*k+4],self.truncation_ratio/c,ax,A,B)

            n,p = len(A)-1,len(B)-1
            for k in range(n): # add the media above the substrate
                ax.text(0,(A[k][0]+A[k+1][0])/2.0,A[k][1],clip_on=True)

            if n>-1:
                ax.text(0,(A[n][0]+self.interface_level)/2.0,A[n][1],clip_on=True)
            
            for k in range(p): # add the media below the substrate
                ax.text(0,(B[k][0]+B[k+1][0])/2.0,B[k][1],clip_on=True)

            if p>-1:
                ax.text(0,(B[p][0]+self.interface_level)/2.0,B[p][1],clip_on=True)

            min_pot = formated_potential.min()
            max_pot = formated_potential.max()
            p = formated_potential.shape[0]/2.0

            if equipotential == True:
                for v in np.linspace(min_pot,max_pot,10):
                    contours = measure.find_contours(formated_potential, v)
                    for n, contour in enumerate(contours):
                        len_contour = len(contour[:,0])
                        ax.plot((contour[:, 1]-p)*x_max/p, (contour[:, 0]-p)*y_max/p, linewidth=1.0,color='k')
                                            
            plt.xlim([-x_max,x_max])
            plt.ylim([-y_max,y_max])
            plt.show()      
            return

        elif output == 'potential_map':

            formated_potential = reformat(self.potential,fmt=fmt,deg=deg)
            
            assert (self.points_file == "classic"), "ERROR: Potential map is only calculated for points_file = classic"
            return formated_potential

        elif output == 'data':

            return self.data

        elif output == 'error':
            
            assert (self.points_file == "surface"), "ERROR: The error is only calculated for points_file = surface"
            return self.error
        
def plot_coating(x_radius,y_radius,media,truncation_ratio,ax,A,B):

# Plot arcs for a particle where its composition is different frow the surrounding media and add to 2 lists where to put the text for the media.
    
    theta = np.arccos(truncation_ratio)*180/np.pi
    
    if media[0] != media[2]:
        coating = Arc((0,0),2.0*x_radius,2.0*y_radius,color='k',linestyle='dashed',linewidth=1.0,theta1=theta+90.0,theta2=90-theta) # add the upper coating
        ax.add_artist(coating)
        A.append([-y_radius,media[2]]) # add the y_extremity of the coating and the material in it.

    if media[1] != media[3]:
        coating = Arc((0,0),2.0*x_radius,2.0*y_radius,color='k',linestyle='dashed',linewidth=1.0,theta2=theta+90.0,theta1=90-theta) # add the lower coating
        ax.add_artist(coating)
        B.append([y_radius,media[3]])        
        
    return
