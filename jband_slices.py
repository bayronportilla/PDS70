import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy import units
from astropy.io import fits
from photutils import EllipticalAnnulus,CircularAnnulus,EllipticalAperture
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_image_jband
import sys
plt.style.use('fancy')

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

def shift(M,c,d):
    a=M.min()
    b=M.max()
    c=-1
    d=+1
    for i in range(0,M.shape[0]):
        for j in range(0,M.shape[1]):
            M[i][j]=c+((d-c)/(b-a))*(M[i][j]-a)
    return M
    

def jband():
    lim=120.0

    ############################################################
    # Absolute paths to files
    run="run_156"
    fits_dir='/data/users/bportilla/runs/final_runs/'+run
    path_fits_image=fits_dir+'/output/RToutObs0001_000001.25.fits.gz'
    path_image_file=fits_dir+'/Image_jband.out'
    path_input_file=fits_dir+'/input.dat'

    ############################################################
    # Fetching information
    imfile=open(path_image_file).readlines()
    for line in imfile:
        if line.split('=')[0]=='MCobs:fov':
            fov=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:npix':
            npix=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:phi':
            phi_image=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:theta':
            theta=float(line.split('=')[1].split('!')[0])
        else:
            continue

    infile=open(path_input_file).readlines()
    for line in infile:
        if line.split('=')[0]=='Distance':
            d=float(line.split('=')[1])


    ############################################################
    # Derived quantities
    pxsize=fov/npix # pixel scale (arcsec/px)
    phi=(phi_image*units.deg).to(units.rad).value # PA from north to east (rad)
    theta=(theta*units.deg).to(units.rad).value # Inclination (rad)
    d=(d*units.pc).to(units.au).value # Distance (au)
    e=np.sin(theta) # eccentricity of the annulus

    ############################################################
    # Load MCMax3D image
    hdulist=fits.open(path_fits_image)
    data_I=hdulist[0].data[0]
    data_Q=hdulist[0].data[1]
    data_U=hdulist[0].data[2]
    data_PI=hdulist[0].data[3]
    
    ############################################################
    # Start here
    #Qphi_g, Uphi_g= combine_Polarizations(data_Q,data_U,((phi)*units.rad).to(units.deg).value)
    phi0=0
    Qphi_g, Uphi_g= combine_Polarizations(data_Q,data_U,phi0)
    #Qphi_g, Uphi_g= combine_Polarizations(data_Q,data_U,0)
    
    # Shift data
    Qphi_g=shift(Qphi_g,-1,1)    
    data_mod=Qphi_g

    ############################################################
    #
    # Get radial profiles - model
    #
    ############################################################
    
    
    ############################################################
    # Input params
    angle_annulus=0#np.pi-phi

    # Determining limit for radial profile
    linear_lim=2*(lim) # AU
    angular_lim=linear_lim/d # rad
    angular_lim=(angular_lim*units.rad).to(units.arcsec).value # arcsec
    pixel_lim=int(round(angular_lim/pxsize))
    xc=0.5*data_mod.shape[0] # Image center in data coordinates
    yc=0.5*data_mod.shape[1] # Image center in data coordinates
    dr=1.0 # Width of the annulus
    a_in_array=[]
    for i in np.arange(yc+dr,yc+0.5*pixel_lim,dr):
        a_in_array.append(i-xc)
    a_out_array=[i+dr for i in a_in_array]
    b_out_array=[i*(1-e**2)**0.5 for i in a_out_array]

    apertures=[EllipticalAnnulus((yc,xc),a_in=ain,a_out=aout,b_out=bout,theta=angle_annulus)
               for (ain,aout,bout) in zip(a_in_array,a_out_array,b_out_array)]

    
    ############################################################
    # Do a check
    a=0.4
    vmin_jband=np.percentile(data_mod,a)
    vmax_jband=np.percentile(data_mod,100-a)
    fig=plt.figure(figsize=(7,6))
    gs=gridspec.GridSpec(2,2,hspace=0.3)
    ax1=plt.subplot(gs[0,0])
    ax2=plt.subplot(gs[0,1])
    ax3=plt.subplot(gs[1,0])
    ax4=plt.subplot(gs[1,1])
    pos=ax1.imshow(data_I,clim=(vmin_jband,vmax_jband))
    ax2.imshow(data_Q,clim=(vmin_jband,vmax_jband))
    ax3.imshow(data_U,clim=(vmin_jband,vmax_jband))
    ax4.imshow(data_PI,clim=(vmin_jband,vmax_jband))
    ax1.set_title(r"$I$")
    ax2.set_title(r"$Q$")
    ax3.set_title(r"$U$")
    ax4.set_title(r"$PI$")
    fig.suptitle(r"$\phi=%.1f$, $\theta=%.1f$, $\phi_0=%.1f$"%(phi_image,theta,phi0))
    plt.show()
    
    
    
    

    return data_mod

jband()
