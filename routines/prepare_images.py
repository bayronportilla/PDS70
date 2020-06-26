from scipy import misc
import scipy
import numpy as np
import sys
from astropy import units
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
from mcmax3d_analysis.mcmax3d_observables import convert_flux,convert_flux_data
from mcmax3d_analysis.mcmax3d_image import display_image
from mcmax3d_analysis.mcmax3d_convolution import convolve_model
from astropy.convolution import Gaussian2DKernel

def combine_Polarizations(Q,U,phi0):
    Qphi= np.zeros(np.shape(Q))
    Uphi= np.zeros(np.shape(Q))
    x0= np.shape(Q)[1]/2#-0.5
    y0= np.shape(Q)[0]/2#-0.5
    phi0=(phi0*units.deg).to(units.rad).value
    for j in range(np.shape(Q)[0]): # over rows
            for i in range(np.shape(Q)[1]): # over columns
                phi= np.arctan2((float(i)-x0),(float(j)-y0))+phi0
                Qphi[i,j]= Q[i,j]*np.cos(2*phi)+U[i,j]*np.sin(2*phi)
                Uphi[i,j]= -Q[i,j]*np.sin(2*phi)+U[i,j]*np.cos(2*phi)
    return Qphi, Uphi


def pivot(data,r,PA,pxsize,d):

    ############################################################
    # Inputs parameters
    xc=data.shape[1]*0.5
    yc=data.shape[0]*0.5
    
    r=(r*units.au).to(units.pc).value # pc
    r=((r/d)*units.rad).to(units.arcsec).value # arcsec
    r=r/pxsize #px
    PA=((PA+90)*units.deg).to(units.rad).value # w.r.t. x-axis 

    xp=r*np.cos(PA)
    yp=r*np.sin(PA)

    x=xp+xc
    y=yp+yc

    return (x,y,data[int(round(y)),int(round(x))])
    

def prepare_model(directory,**kwargs):
    path_input_file=directory+'/input.dat'
    infile=open(path_input_file).readlines()
    for line in infile:
        if line.split('=')[0]=='Distance':
            d=float(line.split('=')[1])

    if kwargs["type"]=="Qphi":
        path_fits_file=directory+"/output/RToutObs0001_000001.25.fits.gz"
        path_image_file=directory+'/Image_jband.out'
            
        ############################################################
        # Loading Image.out info
        imfile=open(path_image_file).readlines()
        for line in imfile:
            if line.split('=')[0]=='MCobs:fov':
                fov=float(line.split('=')[1].split('!')[0])
            elif line.split('=')[0]=='MCobs:npix':
                npix=float(line.split('=')[1].split('!')[0])
            elif line.split('=')[0]=='MCobs:phi':
                phi=float(line.split('=')[1].split('!')[0])
            elif line.split('=')[0]=='MCobs:theta':
                theta=float(line.split('=')[1].split('!')[0])
            else:
                continue

        
        ############################################################
        # Derived quantities
        pxsize=fov/npix # pixel scale (arcsec/px)
        phi=(phi*units.deg).to(units.rad).value # PA from north to east (rad)
        e=np.sin(theta) # eccentricity of the annulus


        ############################################################
        # Load MCMax3D image
        hdulist=fits.open(path_fits_file)
        data_Q=hdulist[0].data[1]
        data_U=hdulist[0].data[2]

    
        ############################################################
        # Creating Qphi, Uphi
        Qphi_g, Uphi_g= combine_Polarizations(data_Q,data_U,0)
        data_mod=Qphi_g
        data_rot=scipy.ndimage.rotate(data_mod,-(158.6-90),reshape=False)
        data=data_rot
        

        ############################################################
        # Searching for flux at a particular pixel 
        r_max_obs=54.64 # AU
        PA_max_obs=158.6 # deg
        xmax=pivot(data,r_max_obs,PA_max_obs,pxsize,d)[0]
        ymax=pivot(data,r_max_obs,PA_max_obs,pxsize,d)[1]
        Bmax=pivot(data,r_max_obs,PA_max_obs,pxsize,d)[2]


        ############################################################
        # Creating info file
        f=open("info_max_Qphi_mod.dat","w")
        f.write("x_max=%d (px)\n"%(xmax))
        f.write("y_max=%d (px)\n"%(ymax))
        f.write("pxsize_mod=%.3f (arcsec/px)\n"%(pxsize))
        f.write("B_max=%.15f (a.u.)\n"%(Bmax))


        ############################################################
        # Creating fits file
        hdu=fits.PrimaryHDU(data)
        hdu.writeto("Qphi_model_rotated.fits",overwrite=True)


    if kwargs["type"]=="alma":
        path_fits_file=directory+"/output/RTout0001_000854.89.fits.gz"
        path_image_file=directory+'/Image_alma.out'
        
        ############################################################
        # Loading Image.out info
        imfile=open(path_image_file).readlines()
        for line in imfile:
            if line.split('=')[0]=='MCobs:fov':
                fov=float(line.split('=')[1].split('!')[0])
            elif line.split('=')[0]=='MCobs:npix':
                npix=float(line.split('=')[1].split('!')[0])
            elif line.split('=')[0]=='MCobs:phi':
                phi=float(line.split('=')[1].split('!')[0])
            elif line.split('=')[0]=='MCobs:theta':
                theta=float(line.split('=')[1].split('!')[0])
            else:
                continue


        ############################################################
        # Derived quantities
        pxsize=fov/npix # pixel scale (arcsec/px)
        phi=(phi*units.deg).to(units.rad).value # PA from north to east (rad)
        e=np.sin(theta) # eccentricity of the annulus


        ############################################################
        # Convolution
        beam_x=0.074 # arcsec
        beam_y=0.057 # arcsec
        beam_angle=63.0 # deg
        hdulist=fits.open(path_fits_file)
        data=hdulist[0].data[0] # mJy/arcsec^2
        data=convolve_model(data,fov,npix,beam_x,beam_y,beam_angle) # mJy/arcsec^2


        ############################################################
        # Convertion from mJy/arcsec^2 to mJy/beam
        beam_area=np.pi*(beam_x)*(beam_y)/(4*np.log(2)) 
        data=data*beam_area # mJy/beam


        ############################################################
        # Rotate image        
        data_rot=scipy.ndimage.rotate(data,-(158.6-90),reshape=False)
        data=data_rot


        ############################################################
        # Searching for flux at a particular pixel 
        r_max_obs=65.44 # AU
        PA_max_obs=326.31 # deg
        xmax=pivot(data,r_max_obs,PA_max_obs,pxsize,d)[0]
        ymax=pivot(data,r_max_obs,PA_max_obs,pxsize,d)[1]
        Bmax=pivot(data,r_max_obs,PA_max_obs,pxsize,d)[2]
        

        ############################################################
        # Creating info file
        f=open("info_max_alma_mod.dat","w")
        f.write("x_max=%d (px)\n"%(xmax))
        f.write("y_max=%d (px)\n"%(ymax))
        f.write("pxsize_mod=%.3f (arcsec/px)\n"%(pxsize))
        f.write("B_max=%.15f (mJy/beam)\n"%(Bmax))
        
        
        # Creating simple file
        ff=open("Bmax.dat","w")
        ff.write("%.15f"%Bmax)
        ff.close()
        # Creating fits file
        hdu=fits.PrimaryHDU(data)
        hdu.writeto("alma_model_rotated.fits",overwrite=True)


    return None


#prepare_model(".",type="alma")
#prepare_model(".",type="Qphi")




