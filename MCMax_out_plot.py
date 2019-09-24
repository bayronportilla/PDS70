import numpy as np
import sys

c=300.0e6 # speed of light m/s

def out_file_MCMax(path_to_file,output_file_name):
    ############################################################
    # MCMax3D data
    data=np.loadtxt(path_to_file)
    x_array=[]
    y_array=[]
    for i in range(0,data.shape[0]):
        factor=c/((data[i][0]*1e-6)**2) # (m/s/m^2)
        flux_min=factor*(data[i][1]*1e-26) # from Jy to W/m^3
        x_array.append(data[i][0])
        y_array.append(flux_min/1e6)
    x_array=np.array(x_array) # microns
    y_array=np.array(y_array) # W/m^2/microns 
    file=open('%s'%output_file_name,'w')
    for i in range(0,len(x_array)):
        file.write('%.16e %.16e\n'%(x_array[i],y_array[i]))
    return file
    
path='/Users/users/bportilla/Documents/first_project/scripts/PDS70/star0001.dat'
out_file_MCMax(path,'converted_star0001.dat')



