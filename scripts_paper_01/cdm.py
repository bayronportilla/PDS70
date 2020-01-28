# Compare Data and Model

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import CubicSpline

def compare_files(file_obs,file_mod):
    ############################################################
    # Import data
    data_obs=np.loadtxt(file_obs)
    data_mod=np.loadtxt(file_mod)

    ############################################################
    # Creating independent arrays
    r_obs=np.reshape(data_obs[:,0:1],data_obs.shape[0])
    f_obs=np.reshape(data_obs[:,1:2],data_obs.shape[0])
    r_mod=np.reshape(data_mod[:,0:1],data_mod.shape[0])
    f_mod=np.reshape(data_mod[:,1:2],data_mod.shape[0])

    ############################################################
    # Prepare cubic spline interpolation
    cs_mod=CubicSpline(r_mod,f_mod)

    ############################################################
    # Interpolate model
    model_interpolated=[]
    radii_interpolated=[]
    for i in np.arange(0,len(r_obs)):
        r=r_obs[i]
        if min(r_mod)<=r<=max(r_mod):
            model_interpolated.append(cs_mod(r))
            radii_interpolated.append(r)

    model_interpolated=np.array(model_interpolated)
    radii_interpolated=np.array(radii_interpolated)

    ############################################################
    # Computing errors
    earray=[abs(i-j)/i*100 for (i,j) in zip(f_obs,model_interpolated)]

    return int(np.ceil(max(earray)))

print(compare_files('../alma_radial_profile_observed.dat','../alma_radial_profile_modeled.dat'))

    
