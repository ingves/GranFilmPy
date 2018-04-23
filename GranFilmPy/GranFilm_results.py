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
    gf.info        =  " Time of generation : %s " % GranFilm_tools.timestamp() 
    
    # Set energy and wavelength
    gf.energy     = data[:,0]
    gf.wavelength = data[:,19]

    # Epsilon
    gf.Epsilon = GranFilm_Results_Epsilon(data_Epsilon,gf)

    # Suceptibilities
    gf.Susceptibilities =  GranFilm_Results_Suceptibilities(data)

    # Polarizabilities
    gf.Polarizabilities =  GranFilm_Results_Polarizabilities(data)
    
    # Reflectivity
    gf.Reflectivity2 = GranFilm_Results_Reflectivity(gf,data_RefTrans)
    
    # Transmittance
    gf.Transmittance2 = GranFilm_Results_Transmittance(gf,data_RefTrans)

    # Potential
    if data_Potential is not None:
        gf.Potential = GranFilm_Results_Potential(gf,data_Potential)
        
        
# --------------- Reflectivity --------------- #

    
class GranFilm_Results_Reflectivity:

    def __init__(self,gf,data):
        
        imu = 1j
        
        # Useful quantities
        medium_1      = gf.param["media"][0]
        medium_2      = gf.param["media"][1]
        n1,n2         = np.sqrt(gf.Epsilon[medium_1]),np.sqrt(gf.Epsilon[medium_2])
        theta_1       = gf.param["theta"]*np.pi/180.0
                             
        # Reflection amplitudes
        self.rp       = data[:,1]  + imu * data[:,2]
        self.rp_flat  = data[:,3]  + imu * data[:,4]
        self.rs       = data[:,5]  + imu * data[:,6]
        self.rs_flat  = data[:,7]  + imu * data[:,8]
        self.rps      = np.sqrt(()
        
        # Reflection coefficients for intensity
        self.Rp       = np.abs(self.rp)**2.0
        self.Rp_flat  = np.abs(self.rp_flat)**2.0
        self.Rs       = np.abs(self.rs)**2.0
        self.Rs_flat  = np.abs(self.rs_flat)**2.0
        self.Rps      = (np.abs(self.rp)**2.0 + np.abs(self.rs)**2.0))/2
        self.Rps_flat = (np.abs(self.rp)**2.0 + np.abs(self.rs)**2.0)/2
        
        # Relative reflection coefficients for intensity
        if gf.param["media"][0] != gf.param["media"][1]:
            self.dRp  = (self.Rp - self.Rp_flat)/self.Rp_flat
            self.dRs  = (self.Rs - self.Rs_flat)/self.Rs_flat
            self.dRps = (self.Rps - self.Rps_flat)/self.Rps_flat
   
    def __call__(self, pol, var='R'):
        
        if var == 'R0':
            s = 'R'+pol+'_flat'
        elif var == 'dR_R':
            s = 'dR'+pol
        else:
            s = var+pol
            
        assert (hasattr(self,s)), "ERROR : var={} and pol={} are not valid arguments".format(var,pol)
        return getattr(self,s)


# --------------- Transmittance --------------- #

    
class GranFilm_Results_Transmittance:

    def __init__(self,gf,data):
        
        # Imaginary unit
        imu = 1j 
                              
        # Tranmission amplitudes
        self.tp       = data[:,9]  + imu * data[:,10]
        self.tp_flat  = data[:,11] + imu * data[:,12]
        self.ts       = data[:,13] + imu * data[:,14]
        self.ts_flat  = data[:,15] + imu * data[:,16]
        
        # Transmission coefficients for intensity
        self.Tp       = np.abs(n2/n1)*(c2/c1)*np.abs(self.tp)**2.0
        self.Tp_flat  = np.abs(n2/n1)*(c2/c1)*np.abs(self.tp_flat)**2.0
        self.Ts       = np.abs(n2/n1)*(c2/c1)*np.abs(self.ts)**2.0
        self.Ts_flat  = np.abs(n2/n1)*(c2/c1)*np.abs(self.ts_flat)**2.0
        self.Tps      = (self.Tp + self.Ts)/2.0
        self.Tps_flat = (self.Tp_flat + self.Ts_flat)/2.0
        
        # Relative transmission coefficients for intensity
        if gf.param["media"][0] != gf.param["media"][1]:
            self.dTp  = (self.Tp - self.Tp_flat)/self.Tp_flat
            self.dTs  = (self.Ts - self.Ts_flat)/self.Ts_flat
            self.dTps = (self.Tps - self.Tps_flat)/self.Tps_flat
   
    def __call__(self, pol, var='T'):
        
        if var == 'T0':
            s = 'T'+pol+'_flat'
        elif var == 'dT_T':
            s = 'dR'+pol
        else:
            s = var+pol
            
        assert (hasattr(self,s)), "ERROR : var={} and pol={} are not valid arguments".format(var,pol)
        return getattr(self,s)


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
                
        imu = 1j                         # imaginary unit
        #
        self.energy = data[:,0]
        #
        self.gamma   = data[:,11]  + imu * data[:,12]
        self.beta    = data[:,13]  + imu * data[:,14]
        #
        self.delta   = data[:,15]  + imu * data[:,16]
        self.tau     = data[:,17]  + imu * data[:,18]


    def __call__(self, type="None", fmt='Complex', deg='True'):
        """
        Get the susceptibility of the system

        Options
        -------
        type: string
            which surface susceptibility to return (gamma, beta, delta or tau)
            Default : None
        fmt : string
           Format of the returned value {'Complex', 'Real', 'Imag', 'Amplitude', 'Phase' }
           Default : Complex
        deg: string
           If fmt=Phase use deg (or not) 
        """           

        # error checking
        FMT={'Complex', 'Real', 'Imag', 'Amplitude', 'Phase' }
        assert (type in {'gamma','beta', 'delta', 'tau'}), "ERROR : type=%s is not supported (valid options: 'gamma','beta','delta','tau') " % type        
        assert (fmt in FMT), "ERROR : type=%s is not supported " % fmt

        if (type=='gamma'):
            susceptibility = self.gamma
        elif (type=='beta'):
            susceptibility = self.beta
        elif (type=='delta'):
            susceptibility = self.delta
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
        self.alpha_10_parallel       = data[:,7]   + imu * data[:,8]
        self.alpha_10_perpendicular  = data[:,9]   + imu * data[:,10]



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
        for material in gf.param["media"]:
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
        self.points_file = gf.param.Potential["points_file"]
        
        if self.points_file == "classic":
            
            x_border = self.data[0][0] # -area_ratio * x_radius
            y_border = self.data[0][2] # -area_ratio * y_radius
            self.x_radius = x_border/-gf.param.Potential["area_ratio_pot"]
            self.y_radius = y_border/-gf.param.Potential["area_ratio_pot"]
            self.area = (x_border,-x_border,-y_border,y_border)
            self.truncation_ratio = gf.param["truncation_ratio"]
            self.interface_level = self.truncation_ratio*self.y_radius
            self.media = gf.param["media"]
            self.radius_ratios = gf.param["radius_ratios"]

            # Reshapment of the data to create the "map" with the values of the potential
            
            n=len(self.data)
            p=int(np.sqrt(n))
            P=np.zeros((p,p),dtype=complex)
            
            for i in range(p):
                for k in range(p):
                    P[k,i]=self.data[i*p+k][4]+1j*self.data[i*p+k][5]
                    
            self.potential = P
            
            x = np.linspace(x_border,-x_border,p)
            y = np.linspace(y_border,-y_border,p)
            self.X,self.Y = np.meshgrid(x,y)

        elif self.points_file == "surface":
            
            # Calculation of the error
            s = 0
            n = len(self.data)
            for k in range(n//2):
                s += np.abs(self.data[2*k][4]-self.data[2*k+1][4])

            self.error = 2*s/n
            
        
    def __call__(self,output='plot',fmt="Real",deg='True',equipotential=True,lines=10,linewidths=1):

        
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

            p = formated_potential.shape[0]/2.0

            if equipotential == True:
                cs = plt.contour(self.X,self.Y,formated_potential,lines,linewidths=linewidths)
                plt.clabel(cs,fontsize=0,inline=0)
                                            
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
