#SJL 5/17
#script to calculate various properties of the surface of a HERCULES planet
#SJL 2/5/18: Corrected and expanded (see notes)

###########################################################
###########################################################
###########################################################
import numpy as np
import sys

#functions for integration
from scipy import integrate



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
    sys.path.insert(0,"C:\\Users\sjl49\Dropbox\Research\Code_repository\Python")


#HERCULES_structures
if sys.platform== 'darwin':
    sys.path.insert(0, path_code+'/Python/SJL_scripts_07_2022')
elif (sys.platform== 'win32') | (sys.platform== 'win64'):
    sys.path.insert(0, 'C:\\Users\sjl49\Dropbox\Research\Code_repository\HERCULES\HERCULES_v1.0\Analysis_scripts')
    
from gradients import *


###############################################################################
###############################################################################
###############################################################################
#Function to calculate the local length per degree and area per degree squared
def calc_local_length_area(p,params,ind_surf):
    #calculated need parameters
    r=p.layers[ind_surf].xi*p.layers[ind_surf].a

    ######
    #lon: is easy
    ######
    dl_lon=r*np.sqrt(1-p.mu**2)

    #######
    #lat: is harder
    #######

    #first calcualte array of angles, indlusing flipped portions to allow continuity at the poles
    theta_combined=np.append(np.append(np.pi-np.arccos(np.flipud(p.mu)),np.arccos(p.mu[1:-1])),-np.arccos(np.flipud(p.mu)))

    #find corresponding r array
    r_combined=np.append(np.append(np.flipud(r),r[1:-1]),np.flipud(r))

    #then find the gradient
    drdtheta=gradient2(theta_combined,r_combined)

    #now find the length but only for the relevant indices
    ind=np.arange(np.size(r), 2*np.size(r))
    dl_lat=np.sqrt(r**2+drdtheta[ind]**2)

    #######
    #Area: is a combination of the two
    #######
    dA=dl_lat*dl_lon

    return dl_lat, dl_lon, dA 

###############################################################################
###############################################################################
###############################################################################
#Function to calculate the distance along the surface between the equator and the pole
def calc_equator_pole_length(p, params, ind_surf):

    #integrate along the surface
    theta=np.arccos(p.layers[ind_surf].mu)
    temp=calc_local_length_area(p,params,ind_surf)
    integrand=temp[0]
    L=-integrate.simps(integrand, theta)

    return L

###############################################################################
###############################################################################
###############################################################################
#Function to calculate the surface area
def calc_surface_area(p, params, ind_surf):

    #integrate the surface
    theta=np.arccos(p.layers[ind_surf].mu)
    temp=calc_local_length_area(p,params,ind_surf)
    integrand=temp[2]
    A=-2.0*2*np.pi*integrate.simps(integrand, theta)

    return A
