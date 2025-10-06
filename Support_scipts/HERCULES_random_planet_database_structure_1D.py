#SJL 3/18
#Modification of HERCULES_random_planet_database_structure_ND to deal with just 1D in L

###############################################################
###############################################################
###############################################################

#import required modules
import numpy as np
import struct
import sys
import os

from copy import copy, deepcopy

#package to use wildcards 
import fnmatch

#interpolator package
from scipy.interpolate import interp1d

if sys.platform== 'darwin':
    
    #get the serial number
    import subprocess
    cmd = "system_profiler SPHardwareDataType | awk '/Serial Number/ {print $4}'"
    result = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True, check=True)
    serial_number = result.stdout.strip().decode('utf-8')
    
    #if uncle bulgaria:
    if serial_number=='D25M61Y9F8JC':
        path_db='/Users/simonlock/Dropbox/'
        path_code='/Users/simonlock/Documents/Bristol_planetary_code_repository'
        path_ody='/Volumes/Lock_onsite_active_memory_extension/Odyssey_transfer_31_3_2019'
    #or if either of the Bristol computers
    elif (serial_number=='C02H91W8PN7C'): #iMac
        path_db='/Users/vq21447/Dropbox/'
        path_code='/Users/vq21447/Documents/Bristol_planetary_code_repository'
        path_ody='/Volumes/Lock_office_active_backup_and_memext/Odyssey_transfer_31_3_2019'
    elif(serial_number=='HNV6D7KP7L'):
        print("CHECK THE LOCATION OF ODYSSEY BACKUP")
        path_db='/Users/vq21447/Dropbox/'
        path_code='/Users/vq21447/Documents/Bristol_planetary_code_repository'
        path_ody='/Volumes/Lock_office_active_backup_and_memext/Odyssey_transfer_31_3_2019'

elif (sys.platform== 'win32') | (sys.platform== 'win64'):
    sys.path.insert(0, 'C:\\Users\sjl49\Dropbox\Research\Code_repository\HERCULES\HERCULES_v1.0\Analysis_scripts')

#HERCULES_structures
if sys.platform== 'darwin':
    sys.path.insert(0, path_code+'/HERCULES/SJL_analysis_scripts_06_2022')
elif (sys.platform== 'win32') | (sys.platform== 'win64'):
    sys.path.insert(0, 'C:\\Users\sjl49\Dropbox\Research\Code_repository\HERCULES\HERCULES_v1.0\Analysis_scripts')
from HERCULES_structures import *
from Mfrac_pressure_calc import *
from surface_size_calc import *
from energy_calculation import *
from calc_gravitational_acceleration import *

if sys.platform== 'darwin':
    sys.path.insert(0, path_code+'/EOS/SJL_scripts_07_2022')
elif (sys.platform== 'win32') | (sys.platform== 'win64'):
    sys.path.insert(0, 'C:\\Users\sjl49\Dropbox\Research\Code_repository\HERCULES\HERCULES_v1.0\Analysis_scripts')
from Gadget_EOS_structure2 import *


if sys.platform== 'darwin':
    sys.path.insert(0, path_code+'/Python/SJL_scripts_07_2022')
elif (sys.platform== 'win32') | (sys.platform== 'win64'):
    sys.path.insert(0, 'C:\\Users\sjl49\Dropbox\Research\Code_repository\HERCULES\HERCULES_v1.0\Analysis_scripts')
from gradients import *

#TESTING
import time


###############################################################
###############################################################
###############################################################
#Structure for a random array of planets with varrying M, L and core fraction
class HERCULES_random_planet_database_1D(object):
    ###################################################
    ###################################################
    #takes as input
    #dir_root: An array of directories that contain various HERCULES run directories. dir_root -> dir_run -> Output
    #run_names: An array of the root names corresponding to the given root directories
    def __init__(self):
        #reference values
        self.ref_L=1.0
        self.min_L=1.0

        #Input bodies
        self.dir_root=[]
        self.run_names=[]

        ## of planets
        self.Np=0
        
        #array of planets
        self.parr=[]
        self.paramsarr=[]

        #Arrays of parameters (e.g., M, fcore etc.)
        self.L=np.empty(0)

        ######CoRoL finding
        #Arrays of parameters + AM for CoRoL and extropolated CoRoL
        self.L_CoRoL=0.0
        self.L_CoRoL_minus1=0.0
        self.L_CoRoL_extrap=0.0

        #Array of the indices of the last and second to last planets
        self.ind_CoRoL=0
        self.ind_CoRoL_minus1=0

        ######Interpolation
        #array of interpolation hulls and information
        self.flag_interp=[]
        self.interp_params=[]
        self.interphull=[]

        #values at the CoRoL
        self.values_CoRoL_extrap=[]
        

    ####################################################
    ####################################################
    ####################################################
    #function to read in an array of planets
    #pass the array of directories and run names for the bodies that need to be included
    #assumes that the max L in each directory marks the limit of the data
    #SJL 8/3/18: Added optional parameters for the dimensions to use and associated parameters (e.g., list of thermal states)
    def make_array(self, dir_root, run_names):
        self.dir_root=deepcopy(dir_root)
        self.run_names=deepcopy(run_names)

        ############################################################
        #extract all sub-CoRoL planets
        #############################################################
        print(self.run_names)
        #check to see if the directory exists
        if os.path.isdir(self.dir_root)==False:
            print(self.run_names)
            print('\t Directory does not exist')
            print('\t '+self.dir_root)
            print('EXITING')
            sys.exit()

        #check to see what successful runs there have been
        runall=os.listdir(self.dir_root+'/')
        runlist=[k for k in runall if fnmatch.fnmatch(k,self.run_names+'*')]
        #print runlist

        #if no runs then let me know
        if np.size(runlist)==0:
            print('\t No matching runs in directory')
            print('EXITING')
            sys.exit()

        #initialise arrays to keep track of the sucessful runs in a directory
        pdir=[]
        paramdir=[]
        Ldir=[]

        #loop over all the runs and find whether they have been successful, failed or have hit the CoRoL
        for k in np.arange(np.size(runlist)):
            dir=self.dir_root+'/'+runlist[k]

            #find the final output for each AM 
            outall=os.listdir(dir+'/Output')
            outlist=[x for x in outall if  fnmatch.fnmatch(x, '*final') ]

            #Only look at successful runs
            if np.size(outlist)!=0:
                #read in final planet
                out=outlist[0]
                #print out
                Hparams=HERCULES_parameters()
                p=HERCULES_planet()
                file = open(dir+'/Output/'+out, "rb")
                Hparams.read_binary(file)
                p.read_binary(file)

                #see if the run was above or below the CoRoL
                if ((((np.absolute(p.Mtot-p.Mtot_tar)/p.Mtot_tar)>1E-13)\
                         |((np.absolute(p.Ltot-p.Ltot_tar)/(p.Ltot_tar+1))>1E-13))\
                        |((((np.amin(p.press)-p.pmin)/p.pmin)>1E-12)\
                              |(np.absolute((p.layers[0].omega-p.omega_rot)/(1+p.omega_rot))>=1E-12))): #not converged properly or above the CoRoL
                    print('\t '+runlist[k])
                    Hparams_fake=HERCULES_parameters()
                    p_fake=p
                    p_fake.Nlayer=np.nan
                    #self.parr.append(p_fake)
                    #self.paramsarr.append(Hparams_fake)
                else: #normal exit, corotating body
                    pdir.append(p)
                    paramdir.append(Hparams)
                    Ldir.append(p.Ltot)
        

        #take the max L value as the limit to the data
        Ldir=np.asarray(Ldir)
        temp=np.argmax(Ldir)
        self.L_CoRoL=Ldir[temp]

        #order the values by AM
        ind=np.argsort(Ldir)
        Ldir=Ldir[ind]

        #add the planets from this directory to the list  
        for mm in ind:
            self.parr.append(pdir[mm])
            self.paramsarr.append(paramdir[mm])


        ########extrapolate to find the CoRoL AM
        #find the indices of the last two AM
        ind=np.argsort(Ldir)
        ind_parr=np.asarray([len(self.parr)-np.size(Ldir)+ind[-2],len(self.parr)-np.size(Ldir)+ind[-1]])
        self.ind_CoRoL=ind_parr[1]
        self.ind_CoRoL_minus1=ind_parr[0]
        self.L_CoRoL_minus1=self.parr[ind_parr[0]].Ltot

        omega_Kep=np.empty(2)
        omega_rot=np.empty(2)
        for mm in np.arange(2):
            #find the radius of all the layers
            radius=np.empty(self.parr[ind_parr[mm]].Nlayer)
            for pp in np.arange(self.parr[ind_parr[mm]].Nlayer):
                radius[pp]=self.parr[ind_parr[mm]].layers[pp].a

            dVdr=gradient2(radius, self.parr[ind_parr[mm]].Ulayers-(((radius*self.parr[ind_parr[mm]].omega_rot)**2)/2.0))
            temp=np.sqrt(-1.0*dVdr/radius)

            omega_Kep[mm]=temp[0]
            omega_rot[mm]=self.parr[ind_parr[mm]].omega_rot

        m_Kep=(omega_Kep[0]-omega_Kep[1])/(Ldir[ind[-2]]-Ldir[ind[-1]])
        m_rot=(omega_rot[0]-omega_rot[1])/(Ldir[ind[-2]]-Ldir[ind[-1]])

        b_Kep=omega_Kep[0]-Ldir[ind[-2]]*m_Kep
        b_rot=omega_rot[0]-Ldir[ind[-2]]*m_rot

        #interp AM at CoRoL
        self.L_CoRoL_extrap=Ldir[ind[-2]]+(Ldir[ind[-2]]-Ldir[ind[-1]])*(omega_Kep[0]-omega_rot[0])/\
            ((omega_rot[0]-omega_rot[1])-(omega_Kep[0]-omega_Kep[1]))
        #omg_CoRoL=m_Kep*AM_CoRoL+b_Kep
            

        ############################################################
        #Now loop and extract the properties of the bodies
        #############################################################
        #define the arrays
        self.Np=np.size(self.parr)
        self.L=np.empty(self.Np)

        #run through the planets and extract the information
        for i in np.arange(self.Np):
            #AM needed in all cases
            self.L[i]=self.parr[i].Ltot

        ############################################################
        #Now normalise to help use
        #############################################################
        #define the reference parameters +AM and normalize axis values and CoRoL values
        #AM
        self.min_L=np.amin(self.L)
        self.ref_L=(np.amax(self.L)-self.min_L)
        self.L=(self.L-self.min_L)/self.ref_L
        self.L_CoRoL=(self.L_CoRoL-self.min_L)/self.ref_L
        self.L_CoRoL_minus1=(self.L_CoRoL_minus1-self.min_L)/self.ref_L
        self.L_CoRoL_extrap=(self.L_CoRoL_extrap-self.min_L)/self.ref_L
        

        return

    ###################################################
    ###################################################
    ###################################################
    #function to initialise the interpolation hulls for the array
    #pass a list of flags for which properties to interpolate and any information required for those
    #0: pcore
    #1: pCMB
    #2: pfmantle (pass the required fraction of mantle. Can be multiple.
    #3: Moment of inertia
    #4: omega
    #5: Kinetic energy #OLD USE FLAG 9
    #6: Grav potential energy #OLD US FLAG 9
    #7: Pressure at a given mass fraction (pass array of mass fractions)
    #8: Density at a given mass fraction (pass array of mass fractions. Assumed same as 7)
    #9: All energy components (Total, Kinetic, Potential, Internal). Requires passing of a Gadget EOS and entropy array for layers. Can also pass an EOS type array, otherwise assumed constant entropy.
    #10: NOT USED
    #11: NOT USED
    #12: The pressure at given eq radius
    #13: The radius of the CMB
    #14: The equatorial and polar radii
    #15: The radius of a given mantle mass fraction
    #16: The local lengths and areas of each mu point
    #17: The surface radius at each mu point
    #18: Effective gravitational acceleration at a given layer
    #######
    #flag_extrap: a flag as to whether to extrapolate beyond data range to interpolated CoRoL
    def initialize_interpolation(self, flag_interp, interp_params, flag_extrap=0):

        #record properties
        self.flag_interp=deepcopy(flag_interp)
        self.interp_params=deepcopy(interp_params)

        #define the size and shape of the interp hull array
        self.interphull=deepcopy(interp_params)

        #set up the CoRoL values
        self.values_CoRoL_extrap=deepcopy(interp_params)

        #array for interpolation (make dimensions depending on dimensions)
        if flag_extrap==0:
            interparr=self.L
            
        elif flag_extrap==1: #if also including CoRoL values add those
            interparr=np.append(self.L,self.L_CoRoL)
            
        else:
            print('WARNING: unidentified flag_extrap')
            print('EXITING')
            sys.exit()


        ###############################################
        #0 core pressure
        ##############################################
        if flag_interp[0]==1:
            p_core=np.empty(self.Np)
            #loop over all the planets
            for i in np.arange(self.Np):
                if np.isnan(self.parr[i].Nlayer):
                    p_core[i]=np.inf#nan
                else:
                    p_core[i]=self.parr[i].pcore

            #if need CoRoL values then loop over all CoRoL points
            if ((flag_extrap==1)):
                m=(p_core[self.ind_CoRoL_minus1]-p_core[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                b=p_core[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                temp=m*self.L_CoRoL_extrap+b

                
                self.values_CoRoL_extrap[0]=temp
                p_core=np.append(p_core,temp)
                
            #set up interpolation hull
            temp=interp1d(interparr,p_core)#, rescale=True)#,fill_value=np.nan)
            self.interphull[0]=temp

        ###############################################
        #1 CMB pressure
        ##############################################
        if flag_interp[1]==1:
            p_CMB=np.empty(self.Np)
            #loop over all the planets
            for i in np.arange(self.Np):
                if np.isnan(self.parr[i].Nlayer):
                    p_CMB[i]=np.nan
                else:
                    core_layer=self.parr[i].Nmaterial-1
                    temp=np.where(self.parr[i].flag_material==core_layer)[0]
                    ind=temp[0]
                    p_CMB[i]=self.parr[i].press[ind]

            
            #if need CoRoL values then loop over all CoRoL points
            if ((flag_extrap==1)):
                m=(p_CMB[self.ind_CoRoL_minus1]-p_CMB[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                b=p_CMB[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                temp=m*self.L_CoRoL_extrap+b

                self.values_CoRoL_extrap[1]=temp
                p_CMB=np.append(p_CMB,temp)
                
            #set up interpolation hull
            temp=interp1d(interparr,p_CMB,fill_value=np.nan)
            self.interphull[1]=temp

        ###############################################
        #2 mantle pressure
        ##############################################
        if flag_interp[2]==1:
            #loop over all the planets
            for i in np.arange(np.size(self.interp_params[2])):
                p_fmantle=np.empty(self.Np)
                for j in np.arange(self.Np):
                    if np.isnan(self.parr[j].Nlayer):
                        p_fmantle[j]=np.nan
                    else:
                        temp=calc_Mfrac_mantle_pressure(self.parr[j], 0.0, self.interp_params[2][i])
                        p_fmantle[j]=temp

                           
                #if need CoRoL values then loop over all CoRoL points
                if ((flag_extrap==1)):
                    m=(p_fmantle[self.ind_CoRoL_minus1]-p_fmantle[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=p_fmantle[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[2][i]=temp
                    p_fmantle=np.append(p_fmantle,temp)


                #set up interpolation hull
                temp=interp1d(interparr,p_fmantle,fill_value=np.nan)
                self.interphull[2][i]=temp

        ###############################################
        #3 Moment of inertia
        ##############################################
        if flag_interp[3]==1:
            I=np.empty(self.Np)
            #loop over all the planets
            for i in np.arange(self.Np):
                if np.isnan(self.parr[i].Nlayer):
                    I[i]=np.nan
                else:
                    I[i]=self.parr[i].I


            #if need CoRoL values then loop over all CoRoL points
            if ((flag_extrap==1)):
                m=(I[self.ind_CoRoL_minus1]-I[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                b=I[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                temp=m*self.L_CoRoL_extrap+b

                self.values_CoRoL_extrap[3]=temp
                I=np.append(I,temp)
                
            #set up interpolation hull
            temp=interp1d(interparr,I,fill_value=np.nan)
            self.interphull[3]=temp

        ###############################################
        #4 Rotation rate
        ##############################################
        if flag_interp[4]==1:
            omg=np.empty(self.Np)
            #loop over all the planets
            for i in np.arange(self.Np):
                if np.isnan(self.parr[i].Nlayer):
                    omg[i]=np.nan
                else:
                    omg[i]=self.parr[i].omega_rot


            #if need CoRoL values then loop over all CoRoL points
            if ((flag_extrap==1)):
                m=(omg[self.ind_CoRoL_minus1]-omg[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                b=omg[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                temp=m*self.L_CoRoL_extrap+b

                self.values_CoRoL_extrap[4]=temp
                omg=np.append(omg,temp)

            #set up interpolation hull
            temp=interp1d(interparr,omg,fill_value=np.nan)
            self.interphull[4]=temp

        ###############################################
        #5+6 Kinetic and potential energy
        ##############################################
        if (flag_interp[5]==1)|(flag_interp[6]==1):
            Ek=np.empty(self.Np)
            Epot=np.empty(self.Np)
            #loop over all the planets
            for i in np.arange(self.Np):
                if np.isnan(self.parr[i].Nlayer):
                    Ek[i]=np.nan
                    Epot[i]=np.nan
                else:
                    temp=calc_pot_kin_energy(self.parr[i], 0)
                    Ek[i]=temp[1]
                    Epot[i]=temp[0]

            #if need CoRoL values then loop over all CoRoL points
            if ((flag_extrap==1)):
                #kinetic
                m=(Ek[self.ind_CoRoL_minus1]-Ek[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                b=Ek[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                temp=m*self.L_CoRoL_extrap+b

                Ek_CoRoL=temp

                #potential
                m=(Epot[self.ind_CoRoL_minus1]-Epot[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                b=Epot[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                temp=m*self.L_CoRoL_extrap+b

                Epot_CoRoL=temp

                self.values_CoRoL_extrap[5]=Ek_CoRoL
                self.values_CoRoL_extrap[6]=Epot_CoRoL

                Ek=np.append(Ek,Ek_CoRoL)
                Epot=np.append(Epot,Epot_CoRoL)
                
            #set up interpolation hull
            if (flag_interp[5]==1):  
                temp=interp1d(interparr,Ek,fill_value=np.nan)
                self.interphull[5]=temp
            if (flag_interp[6]==1):  
                temp=interp1d(interparr,Epot,fill_value=np.nan)
                self.interphull[6]=temp

        ###############################################
        #7+8 Pressure and density at a mass fraction
        ##############################################
        if (flag_interp[7]==1)|(flag_interp[8]==1):
            #loop over all the planets
            for i in np.arange(np.size(self.interp_params[7])):
                p_f=np.empty(self.Np)
                rho_f=np.empty(self.Np)
                for j in np.arange(self.Np):
                    if np.isnan(self.parr[j].Nlayer):
                        p_f[j]=np.nan
                        rho_f[j]=np.nan
                    else:
                        temp=calc_Mfrac_pressure_rho(self.parr[j], 0.0, self.interp_params[7][i])
                        p_f[j]=temp[0]
                        rho_f[j]=temp[1]

                #if need CoRoL values then loop over all CoRoL points
                if ((flag_extrap==1)):
                    #pressure
                    m=(p_f[self.ind_CoRoL_minus1]-p_f[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=p_f[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    p_f_CoRoL=temp

                    #density
                    m=(rho_f[self.ind_CoRoL_minus1]-rho_f[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=rho_f[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    rho_f_CoRoL=temp


                    self.values_CoRoL_extrap[7][i]=p_f_CoRoL
                    self.values_CoRoL_extrap[8][i]=rho_f_CoRoL
                    p_f=np.append(p_f,p_f_CoRoL)
                    rho_f=np.append(rho_f,rho_f_CoRoL)


                #set up interpolation hull
                if (flag_interp[7]==1):
                    temp=interp1d(interparr,p_f,fill_value=np.nan)
                    self.interphull[7][i]=temp
                if (flag_interp[8]==1):
                    temp=interp1d(interparr,rho_f,fill_value=np.nan)
                    self.interphull[8][i]=temp

        
        ###############################################
        #9 All energy components
        ##############################################
        if (len(flag_interp)>=10):
            if (flag_interp[9]==1):
                #extract parameters
                EOS=interp_params[9][0]
                Sarr=interp_params[9][1]
                #if have been given EOS types then use it
                if np.shape(interp_params[9])[0]>2:
                    EOS_type=interp_params[9][2]

                Nmat=np.size(Sarr)

                #reshape the interp hull array
                self.interphull[9]=[0,0,0,[0.0]*Nmat]
                self.values_CoRoL_extrap[9]=[0,0,0,[0.0]*Nmat]

                #initialize arrays
                Etot=np.empty(self.Np)
                Ek=np.empty(self.Np)
                Epot=np.empty(self.Np)
                U=np.empty([self.Np, Nmat])
                #loop over all the planets
                for i in np.arange(self.Np):
                    if np.isnan(self.parr[i].Nlayer):
                        Etot[i]=np.nan
                        Ek[i]=np.nan
                        Epot[i]=np.nan
                        for j in np.arange(Nmat):
                            U[i][j]=np.nan
                    else:
                        if np.shape(interp_params[9])[0]>2:
                            temp=calc_pot_kin_int_energy_EOStype(self.parr[i], 0,EOS,Sarr,EOS_type=EOS_type)
                        else:
                            temp=calc_pot_kin_int_energy_EOStype(self.parr[i], 0,EOS,Sarr)
                        Etot[i]=np.sum(np.hstack(temp))
                        Ek[i]=temp[0]
                        Epot[i]=temp[1]
                        for j in np.arange(Nmat):
                            U[i][j]=temp[2][j]


                #if need CoRoL values then loop over all CoRoL points
                if ((flag_extrap==1)):

                    #Total
                    m=(Etot[self.ind_CoRoL_minus1]-Etot[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=Etot[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    Etot_CoRoL=temp

                    #Kinetic
                    m=(Ek[self.ind_CoRoL_minus1]-Ek[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=Ek[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    Ek_CoRoL=temp

                    #Potential
                    m=(Epot[self.ind_CoRoL_minus1]-Epot[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=Epot[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    Epot_CoRoL=temp

                    self.values_CoRoL_extrap[9][0]=Etot_CoRoL
                    self.values_CoRoL_extrap[9][1]=Ek_CoRoL
                    self.values_CoRoL_extrap[9][2]=Epot_CoRoL

                    Etot=np.append(Etot,Etot_CoRoL)
                    Ek=np.append(Ek,Ek_CoRoL)
                    Epot=np.append(Epot,Epot_CoRoL)

                    #seperate loop for all the different layers for internal energy
                    U_CoRoL=np.empty([1,Nmat])

                    #For each layer 
                    for j in np.arange(Nmat):
                        m=(U[self.ind_CoRoL_minus1][j]-U[self.ind_CoRoL][j])/(self.L_CoRoL_minus1-self.L_CoRoL)
                        b=U[self.ind_CoRoL_minus1][j]-self.L_CoRoL_minus1*m
                        temp=m*self.L_CoRoL_extrap+b

                        U_CoRoL[0,j]=temp

                    for j in np.arange(Nmat):
                        self.values_CoRoL_extrap[9][3][j]=U_CoRoL[0,j]

                    U=np.append(U,U_CoRoL,0)

                #set up interpolation hull
                temp=interp1d(interparr,Etot,fill_value=np.nan)
                self.interphull[9][0]=temp
                temp=interp1d(interparr,Ek,fill_value=np.nan)
                self.interphull[9][1]=temp
                temp=interp1d(interparr,Epot,fill_value=np.nan)
                self.interphull[9][2]=temp
                for j in np.arange(Nmat):
                    temp=interp1d(interparr,U[:,j],fill_value=np.nan)
                    self.interphull[9][3][j]=temp


        ###############################################
        #12 Pressure at given eq and polar radii
        ##############################################
        if (len(flag_interp)>=13):
            if flag_interp[12]==1:
                self.interphull[12]=np.empty([np.size(self.interp_params[12]), 2]).tolist()
                self.values_CoRoL_extrap[12]=np.empty([np.size(self.interp_params[12]), 2]).tolist()
                #loop over all the planets
                for i in np.arange(np.size(self.interp_params[12])):
                    p_Req=np.empty(self.Np)
                    p_Rpl=np.empty(self.Np)
                    for j in np.arange(self.Np):
                        if np.isnan(self.parr[j].Nlayer):
                            p_Req[j]=np.nan
                            p_Rpl[j]=np.nan
                        else:
                            temp=calc_press_rad(self.parr[j], 0.0, self.interp_params[12][i])
                            p_Req[j]=temp[0]
                            temp=calc_press_rad_pole(self.parr[j], 0.0, self.interp_params[12][i])
                            p_Rpl[j]=temp[0]


                    #if need CoRoL values then loop over all CoRoL points
                    if ((flag_extrap==1)):
                        m=(p_Req[self.ind_CoRoL_minus1]-p_Req[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                        b=p_Req[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                        temp=m*self.L_CoRoL_extrap+b

                        self.values_CoRoL_extrap[12][i]=temp
                        p_Req=np.append(p_Req,temp)

                        m=(p_Rpl[self.ind_CoRoL_minus1]-p_Rpl[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                        b=p_Rpl[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                        temp=m*self.L_CoRoL_extrap+b

                        self.values_CoRoL_extrap[12][i]=temp
                        p_Rpl=np.append(p_Rpl,temp)


                    #set up interpolation hull
                    temp=interp1d(interparr,p_Req,fill_value=np.nan)
                    self.interphull[12][i][0]=temp

                    temp=interp1d(interparr,p_Rpl,fill_value=np.nan)
                    self.interphull[12][i][1]=temp


        ###############################################
        #13 Radius of CMB
        ##############################################
        if (len(flag_interp)>=14):
            if flag_interp[13]==1:
                rCMB_eq=np.empty(self.Np)
                rCMB_pl=np.empty(self.Np)
                #loop over all the planets
                for i in np.arange(self.Np):
                    if np.isnan(self.parr[i].Nlayer):
                        rCMB_eq[i]=np.nan
                        rCMB_pl[i]=np.nan
                    else:
                        core_layer=self.parr[i].Nmaterial-1
                        temp=np.where(self.parr[i].flag_material==core_layer)[0]
                        ind=temp[0]
                        rCMB_eq[i]=self.parr[i].layers[ind].a
                        rCMB_pl[i]=self.parr[i].layers[ind].b


                #if need CoRoL values then loop over all CoRoL points
                if ((flag_extrap==1)):
                    m=(rCMB_eq[self.ind_CoRoL_minus1]-rCMB_eq[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=rCMB_eq[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[1]=temp
                    rCMB_eq=np.append(rCMB_eq,temp)

                    m=(rCMB_pl[self.ind_CoRoL_minus1]-rCMB_pl[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=rCMB_pl[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[1]=temp
                    rCMB_pl=np.append(rCMB_pl,temp)

                #set up interpolation hull
                temp=interp1d(interparr,rCMB_eq,fill_value=np.nan)
                self.interphull[13][0]=temp

                temp=interp1d(interparr,rCMB_pl,fill_value=np.nan)
                self.interphull[13][1]=temp

        

        ###############################################
        #14 Polar and equatorial radii
        ##############################################
        if (len(flag_interp)>=15):
            if flag_interp[14]==1:
                amax=np.empty(self.Np)
                bmax=np.empty(self.Np)
                #loop over all the planets
                for i in np.arange(self.Np):
                    if np.isnan(self.parr[i].Nlayer):
                        amax[i]=np.nan
                        bmax[i]=np.nan
                    else:
                        amax[i]=self.parr[i].amax
                        bmax[i]=self.parr[i].amax*self.parr[i].aspect


                #if need CoRoL values then loop over all CoRoL points
                if ((flag_extrap==1)):
                    m=(amax[self.ind_CoRoL_minus1]-amax[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=amax[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[1]=temp
                    amax=np.append(amax,temp)

                    m=(bmax[self.ind_CoRoL_minus1]-bmax[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=bmax[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[1]=temp
                    bmax=np.append(bmax,temp)

                #set up interpolation hull
                temp=interp1d(interparr,amax,fill_value=np.nan)
                self.interphull[14][0]=temp
                
                temp=interp1d(interparr,bmax,fill_value=np.nan)
                self.interphull[14][1]=temp


        
        ###############################################
        #15 mantle mass fraction radii
        ##############################################
        if (len(flag_interp)>=16):
            if flag_interp[15]==1:

                #reshape the interp hull array
                self.interphull[15]=np.empty([np.size(self.interp_params[15]), 2]).tolist()
                self.values_CoRoL_extrap[15]=np.empty([np.size(self.interp_params[15]), 2]).tolist()

                #loop over all the planets
                for i in np.arange(np.size(self.interp_params[15])):
                    a_fmantle=np.empty(self.Np)
                    b_fmantle=np.empty(self.Np)
                    for j in np.arange(self.Np):
                        if np.isnan(self.parr[j].Nlayer):
                            a_fmantle[j]=np.nan
                            b_fmantle[j]=np.nan
                        else:
                            temp=calc_Mfrac_mantle_rad(self.parr[j], 0.0, self.interp_params[15][i])
                            a_fmantle[j]=temp[0]
                            b_fmantle[j]=temp[1]


                    #if need CoRoL values then loop over all CoRoL points
                    if ((flag_extrap==1)):
                        m=(a_fmantle[self.ind_CoRoL_minus1]-a_fmantle[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                        b=a_fmantle[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                        temp=m*self.L_CoRoL_extrap+b

                        self.values_CoRoL_extrap[15][i][0]=temp
                        a_fmantle=np.append(a_fmantle,temp)

                        m=(b_fmantle[self.ind_CoRoL_minus1]-b_fmantle[self.ind_CoRoL])/(self.L_CoRoL_minus1-self.L_CoRoL)
                        b=b_fmantle[self.ind_CoRoL_minus1]-self.L_CoRoL_minus1*m
                        temp=m*self.L_CoRoL_extrap+b

                        self.values_CoRoL_extrap[15][i][1]=temp
                        b_fmantle=np.append(b_fmantle,temp)

                    #set up interpolation hull
                    temp=interp1d(interparr,a_fmantle,fill_value=np.nan)
                    self.interphull[15][i][0]=temp
                    temp=interp1d(interparr,b_fmantle,fill_value=np.nan)
                    self.interphull[15][i][1]=temp

        ###############################################
        #16 Length/area of surface
        ##############################################
        if (len(flag_interp)>=17):
            if flag_interp[16]==1:

                #reshape the interp hull array
                self.interphull[16]=np.empty([self.parr[0].Nmu,6]).tolist()
                self.values_CoRoL_extrap[16]=np.empty([self.parr[0].Nmu,6]).tolist()
                
                llat=np.empty((self.Np,self.parr[0].Nmu))
                llon=np.empty((self.Np,self.parr[0].Nmu))
                A=np.empty((self.Np,self.parr[0].Nmu))
                #loop over all the planets
                for i in np.arange(self.Np):
                    if np.isnan(self.parr[i].Nlayer):
                        llat[i,:]=np.ones(self.parr[0].Nmu)*np.nan
                        llon[i,:]=np.ones(self.parr[0].Nmu)*np.nan
                        A[i,:]=np.ones(self.parr[0].Nmu)*np.nan
                    else:
                        (llat[i,:],llon[i,:],A[i,:])=calc_local_length_area(self.parr[i],self.paramsarr[i],0)

                #also calculate the rate of change with respect to AM
                dllat_dL=np.empty((self.Np,self.parr[0].Nmu))
                dllon_dL=np.empty((self.Np,self.parr[0].Nmu))
                dA_dL=np.empty((self.Np,self.parr[0].Nmu))
                for i in np.arange(self.parr[0].Nmu):
                    dllat_dL[:,i]=gradient2(self.L,llat[:,i])
                    dllon_dL[:,i]=gradient2(self.L,llon[:,i])
                    dA_dL[:,i]=gradient2(self.L,A[:,i])


                #if need CoRoL values then loop over all CoRoL points
                if ((flag_extrap==1)):
                    m=(llat[self.ind_CoRoL_minus1,:]-llat[self.ind_CoRoL,:])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=llat[self.ind_CoRoL_minus1,:]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[16][:][0]=temp
                    llat=np.row_stack((llat,np.reshape(temp,(1,self.parr[0].Nmu))))

                    m=(llon[self.ind_CoRoL_minus1,:]-llon[self.ind_CoRoL,:])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=llon[self.ind_CoRoL_minus1,:]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[16][:][1]=temp
                    llon=np.row_stack((llon,np.reshape(temp,(1,self.parr[0].Nmu))))

                    m=(A[self.ind_CoRoL_minus1,:]-A[self.ind_CoRoL,:])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=A[self.ind_CoRoL_minus1,:]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[16][:][2]=temp
                    A=np.row_stack((A,np.reshape(temp,(1,self.parr[0].Nmu))))

                    #
                    m=(dllat_dL[self.ind_CoRoL_minus1,:]-dllat_dL[self.ind_CoRoL,:])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=dllat_dL[self.ind_CoRoL_minus1,:]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[16][:][3]=temp
                    dllat_dL=np.row_stack((dllat_dL,np.reshape(temp,(1,self.parr[0].Nmu))))

                    m=(dllon_dL[self.ind_CoRoL_minus1,:]-dllon_dL[self.ind_CoRoL,:])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=dllon_dL[self.ind_CoRoL_minus1,:]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[16][:][4]=temp
                    dllon_dL=np.row_stack((dllon_dL,np.reshape(temp,(1,self.parr[0].Nmu))))

                    m=(dA_dL[self.ind_CoRoL_minus1,:]-dA_dL[self.ind_CoRoL,:])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=dA_dL[self.ind_CoRoL_minus1,:]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[16][:][5]=temp
                    dA_dL=np.row_stack((dA_dL,np.reshape(temp,(1,self.parr[0].Nmu))))

                #set up interpolation hull
                for i in np.arange(self.parr[0].Nmu):
                    temp=interp1d(interparr,llat[:,i],fill_value=np.nan)
                    self.interphull[16][i][0]=temp
                
                    temp=interp1d(interparr,llon[:,i],fill_value=np.nan)
                    self.interphull[16][i][1]=temp

                    temp=interp1d(interparr,A[:,i],fill_value=np.nan)
                    self.interphull[16][i][2]=temp

                    #
                    temp=interp1d(interparr,dllat_dL[:,i],fill_value=np.nan)
                    self.interphull[16][i][3]=temp
                
                    temp=interp1d(interparr,dllon_dL[:,i],fill_value=np.nan)
                    self.interphull[16][i][4]=temp

                    temp=interp1d(interparr,dA_dL[:,i],fill_value=np.nan)
                    self.interphull[16][i][5]=temp


        ###############################################
        #17 Radius of surface
        ##############################################
        if (len(flag_interp)>=18):
            if flag_interp[17]==1:

                #reshape the interp hull array
                self.interphull[17]=np.empty(self.parr[0].Nmu).tolist()
                self.values_CoRoL_extrap[17]=np.empty(self.parr[0].Nmu).tolist()
                
                rsurf=np.empty((self.Np,self.parr[0].Nmu))
                #loop over all the planets
                for i in np.arange(self.Np):
                    if np.isnan(self.parr[i].Nlayer):
                        rsurf[i,:]=np.ones(self.parr[i].Nlayer)*np.nan
                    else:
                        rsurf[i,:]=self.parr[i].layers[0].xi*self.parr[i].amax
                        

                #if need CoRoL values then loop over all CoRoL points
                if ((flag_extrap==1)):

                    m=(rsurf[self.ind_CoRoL_minus1,:]-rsurf[self.ind_CoRoL,:])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=rsurf[self.ind_CoRoL_minus1,:]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[17][:]=temp
                    rsurf=np.row_stack((rsurf,np.reshape(temp,(1,self.parr[0].Nmu))))

                #set up interpolation hull
                for i in np.arange(self.parr[0].Nmu):

                    temp=interp1d(interparr,rsurf[:,i],fill_value=np.nan)
                    self.interphull[17][i]=temp


        ###############################################
        #18 Effective gravity at a given layer
        ##############################################
        if (len(flag_interp)>=19):
            if flag_interp[18]==1:

                #reshape the interp hull array
                self.interphull[18]=np.empty([self.parr[0].Nmu,3]).tolist()
                self.values_CoRoL_extrap[18]=np.empty([self.parr[0].Nmu,3]).tolist()
                
                gtot=np.empty((self.Np,self.parr[0].Nmu))
                grxy=np.empty((self.Np,self.parr[0].Nmu))
                gz=np.empty((self.Np,self.parr[0].Nmu))
                #loop over all the planets
                for i in np.arange(self.Np):
                    if np.isnan(self.parr[i].Nlayer):
                        gtot[i,:]=np.ones(self.parr[0].Nmu)*np.nan
                        grxy[i,:]=np.ones(self.parr[0].Nmu)*np.nan
                        gz[i,:]=np.ones(self.parr[0].Nmu)*np.nan
                    else:
                        (grxy[i,:],gz[i,:])=calc_total_effective_gravity_2ndorder(self.parr[i],self.paramsarr[i],self.interp_params[18])

                gtot=np.sqrt(grxy**2+gz**2)

                #if need CoRoL values then loop over all CoRoL points
                if ((flag_extrap==1)):
                    m=(gtot[self.ind_CoRoL_minus1,:]-gtot[self.ind_CoRoL,:])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=gtot[self.ind_CoRoL_minus1,:]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[18][:][0]=temp
                    gtot=np.row_stack((gtot,np.reshape(temp,(1,self.parr[0].Nmu))))

                    m=(grxy[self.ind_CoRoL_minus1,:]-grxy[self.ind_CoRoL,:])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=grxy[self.ind_CoRoL_minus1,:]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[18][:][1]=temp
                    grxy=np.row_stack((grxy,np.reshape(temp,(1,self.parr[0].Nmu))))

                    m=(gz[self.ind_CoRoL_minus1,:]-gz[self.ind_CoRoL,:])/(self.L_CoRoL_minus1-self.L_CoRoL)
                    b=gz[self.ind_CoRoL_minus1,:]-self.L_CoRoL_minus1*m
                    temp=m*self.L_CoRoL_extrap+b

                    self.values_CoRoL_extrap[18][:][2]=temp
                    gz=np.row_stack((gz,np.reshape(temp,(1,self.parr[0].Nmu))))

                #set up interpolation hull
                for i in np.arange(self.parr[0].Nmu):
                    temp=interp1d(interparr,gtot[:,i],fill_value=np.nan)
                    self.interphull[18][i][0]=temp
                
                    temp=interp1d(interparr,grxy[:,i],fill_value=np.nan)
                    self.interphull[18][i][1]=temp

                    temp=interp1d(interparr,gz[:,i],fill_value=np.nan)
                    self.interphull[18][i][2]=temp
        

        return

    ###################################################
    ###################################################
    ###################################################
    #interpolate the desired properties
    #interparr passes an array of values to interpolate but only ones that have been initialized
    #pass a list of flags for which properties to interpolate and any information required for those
    #0: pcore
    #1: pCMB
    #2: pfmantle (pass the required fraction of mantle. Can be multiple.
    #3: Moment of inertia
    #4: omega
    #5: Kinetic energy
    #6: Grav potential energy
    #7: Pressure at a given mass fraction (pass array of mass fractions)
    #8: Density at a given mass fraction (pass array of mass fractions. Assumed same as 7)
    #9: All energy components (Total, Kinetic, Potential, Internal). Requires passing of a Gadget EOS and entropy array for layers
    #10: NOT USED
    #11: NOT USED
    #12: The pressure at given eq and pl radius
    #13: The radius of the CMB
    #14: The equatorial and polar radii
    #15: The radius of a given mantle mass fraction
    #16: Local lenths and area at each mu point
    #17: The surface radius at each mu point
    #18: The effective gravity at a given layer
    #######
    #flag_extrap: a flag as to whether to extrapolate beyond data range to interpolated CoRoL
    def interp_database(self, L, flag_interp, interp_params, flag_extrap=0):
        #initialize output array
        output=deepcopy(interp_params)     

        #convert input to normalized values
        L=(L-self.min_L)/self.ref_L

        ##################################
        ##################################
        #interpolate to find whether the point is above the limit of the data (CoRoL)
        if flag_extrap==0:
            Ltest=self.L_CoRoL
        elif flag_extrap==1:
            Ltest=self.L_CoRoL_extrap
        else:
            print('ERROR: Unknown flag_extrap')
            print('EXITING')
            sys.exit()

        #if above then return nan and exit
        if L>Ltest:
            output[0]=np.nan
            output[1]=np.nan
            for i in np.arange(np.size(interp_params[2])):
                output[2][i]=np.nan
            output[3]=np.nan
            output[4]=np.nan
            output[5]=np.nan
            output[6]=np.nan
            for i in np.arange(np.size(interp_params[7])):
                output[7][i]=np.nan
            for i in np.arange(np.size(interp_params[8])):
                output[8][i]=np.nan

            Sarr=interp_params[9][1]
            Nmat=np.size(Sarr)
            output[9]=[0,0,0,[0.0]*Nmat]
            output[9][0]=np.nan
            output[9][1]=np.nan
            output[9][2]=np.nan
            for i in np.arange(Nmat):
                output[9][3][i]=np.nan

            return output
    
        #interp all properties using the pre-initialized hulls
        if flag_interp[0]==1:
            output[0]=self.interphull[0]((L))
        #
        if flag_interp[1]==1:
            output[1]=self.interphull[1]((L))
        #
        if flag_interp[2]==1:
            for i in np.arange(np.size(interp_params[2])):
                output[2][i]=self.interphull[2][i]((L))
        #
        if flag_interp[3]==1:
            output[3]=self.interphull[3]((L))
        #
        if flag_interp[4]==1:
            output[4]=self.interphull[4]((L))
        #
        if flag_interp[5]==1:
            output[5]=self.interphull[5]((L))
        #
        if flag_interp[6]==1:
            output[6]=self.interphull[6]((L))
        #
        if flag_interp[7]==1:
            for i in np.arange(np.size(interp_params[7])):
                output[7][i]=self.interphull[7][i]((L))
        #
        if flag_interp[8]==1:
            for i in np.arange(np.size(interp_params[7])):
                output[8][i]=self.interphull[8][i]((L))
        #
        if (len(flag_interp)>=10):
            if (flag_interp[9]==1):
                Sarr=interp_params[9][1]
                Nmat=np.size(Sarr)

                output[9]=[0,0,0,[0.0]*Nmat]
                output[9][0]=self.interphull[9][0]((L))
                output[9][1]=self.interphull[9][1]((L))
                output[9][2]=self.interphull[9][2]((L))
                for j in np.arange(Nmat):
                    output[9][3][j]=self.interphull[9][3][j]((L))
        

        #
        if (len(flag_interp)>=13):
            if flag_interp[12]==1:
                output[12]=np.empty([np.size(self.interp_params[12]), 2]).tolist()
                for i in np.arange(np.size(interp_params[12])):
                    output[12][i][0]=self.interphull[12][i][0]((L))
                    output[12][i][1]=self.interphull[12][i][1]((L))
            #
        if (len(flag_interp)>=14):
            if flag_interp[13]==1:
                output[13][0]=self.interphull[13][0]((L))
                output[13][1]=self.interphull[13][1]((L))

        if (len(flag_interp)>=15):
            if flag_interp[14]==1:
                output[14][0]=self.interphull[14][0]((L))
                output[14][1]=self.interphull[14][1]((L))
        #
        if (len(flag_interp)>=16):
            if flag_interp[15]==1:
                output[15]=np.empty([np.size(self.interp_params[15]), 2]).tolist()
                for i in np.arange(np.size(interp_params[15])):
                    output[15][i][0]=self.interphull[15][i][0]((L))
                    output[15][i][1]=self.interphull[15][i][1]((L))


        #
        if (len(flag_interp)>=17):
            if flag_interp[16]==1:
                output[16]=np.empty([self.parr[0].Nmu, 6]).tolist()
                for i in np.arange(self.parr[0].Nmu):
                    output[16][i][0]=self.interphull[16][i][0]((L))
                    output[16][i][1]=self.interphull[16][i][1]((L))
                    output[16][i][2]=self.interphull[16][i][2]((L))
                    output[16][i][3]=self.interphull[16][i][3]((L))/self.ref_L
                    output[16][i][4]=self.interphull[16][i][4]((L))/self.ref_L
                    output[16][i][5]=self.interphull[16][i][5]((L))/self.ref_L

        #
        if (len(flag_interp)>=18):
            if flag_interp[17]==1:
                output[17]=np.empty(self.parr[0].Nmu).tolist()
                for i in np.arange(self.parr[0].Nmu):
                    output[17][i]=self.interphull[17][i]((L))


        #
        if (len(flag_interp)>=19):
            if flag_interp[18]==1:
                output[18]=np.empty([self.parr[0].Nmu, 3]).tolist()
                for i in np.arange(self.parr[0].Nmu):
                    output[18][i][0]=self.interphull[18][i][0]((L))
                    output[18][i][1]=self.interphull[18][i][1]((L))
                    output[18][i][2]=self.interphull[18][i][2]((L))
        
        #sys.exit()

        return output
