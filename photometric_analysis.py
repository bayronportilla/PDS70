import numpy as np
from astropy import constants as cte
import matplotlib.pyplot as plt
import sys

############################################################
# Loading data
photometric_data_keppler=np.loadtxt('PDS70_photometry_keppler.dat',skiprows=2,usecols=(0,1,2))
x_photometry=np.reshape(photometric_data_keppler[:,0:1],photometric_data_keppler.shape[0]) # microns
y_photometry=np.reshape(photometric_data_keppler[:,1:2],photometric_data_keppler.shape[0]) # mJy
y_photometry_error=np.reshape(photometric_data_keppler[:,2:3],photometric_data_keppler.shape[0]) # mJy

############################################################
# Converting mJy to W/m^3
x_array=[]
y_array=[]
y_error_array=[]
for i in range(0,len(x_photometry)):
    factor=cte.c.value/((x_photometry[i]*1e-6)**2) # (m/s/m^2)
    flux=factor*((y_photometry[i]/1e3)*1e-26) # from mJy to W/m^3
    error=factor*((y_photometry_error[i]/1e3)*1e-26) # from mJy to W/m^3
    x_array.append(x_photometry[i]) # microns
    y_array.append(flux/1e6) # W/m^2/microns
    y_error_array.append(error/1e6) # W/m^2/microns

############################################################
# Creating file
file=open('PDS70_photometry.dat','w')
for i in range(0,len(x_array)):
    file.write('%.15e %.15e %.15e\n'%(x_array[i],y_array[i],y_error_array[i]))


