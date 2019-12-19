import numpy as np
import sys
from astropy import units
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
from photutils import EllipticalAnnulus
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_image
from matplotlib.patches import Ellipse
plt.style.use('fancy')

def alma_data(rundir,fits_image,mode,lim):

    ############################################################
    # Absolute paths to files
    path_fits_file=rundir+'output/'+fits_image
    path_image_file=rundir+'Image.out'
    path_input_file=rundir+'input.dat'


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


    ############################################################
    # For photometry
    e=np.sin(theta) # eccentricity of the annulus


    ############################################################
    # Load data
    hdulist=fits.open(path_fits_file)
    assert(mode in ['alma','jband']),"%s This mode is not recognized"%mode
    if mode=='alma':
        data=hdulist[0].data[0]
        data=convolve_image(data)
    elif mode=='jband':
        data_Q=hdulist[0].data[1]
        data_U=hdulist[0].data[2]
        data=data_Q*np.cos(2*phi) + data_U*np.sin(2*phi)


    ############################################################
    # Brightness profile
    linear_lim=2*(lim) # AU
    angular_lim=((linear_lim/d)*units.rad).to(units.arcsec).value # arcsec
    pixel_lim=int(round(angular_lim/pxsize))

    # Creating elliptical aperture object in pixel coordinates
    xc=0.5*data.shape[0]
    yc=0.5*data.shape[1]
    dr=1.0
    ain_array=[]
    for i in np.arange(xc+dr,xc+0.5*pixel_lim,dr):
        ain_array.append(i-xc)
    aout_array=[i+dr for i in ain_array]
    bout_array=[i*(1-e**2)**0.5 for i in aout_array]

    # Creating apertures
    apertures=[EllipticalAnnulus((xc,yc),a_in=ain,a_out=aout,b_out=bout,theta=0.0)
               for (ain,aout,bout) in zip(ain_array,aout_array,bout_array)]
   
    
    ############################################################
    # Do a check?
    a=0.01
    vmin=np.percentile(data,a)
    vmax=np.percentile(data,100-a)
    aperture=apertures[-1]
    fig=plt.figure()
    ax=plt.axes()
    ax.imshow(data,clim=(vmin,vmax))
    aperture.plot(color='red',lw=1)
    if mode=="alma":
        fig.savefig("alma_image.png")
    elif mode=="jband":
        fig.savefig("jband_image.png")
    #plt.show()
    

    # Do photometry
    phot_table=aperture_photometry(data,apertures)
    col_values=[]
    for col in phot_table.colnames:
        col_values.append(phot_table[col][0])
    brightness=[col_values[i] for i in range(3,len(col_values))]
    brightness=np.array(brightness)

    #if mode=='jband':

    # Create arrays for output file
    r_arcsec=[(j+0.5*(i-j))*pxsize for (i,j) in zip(aout_array,ain_array)] # arcsec
    r_rad=[(i*units.arcsec).to(units.rad).value for i in r_arcsec] # rad
    r_au=[((i*d)*units.au).value for i in r_rad] # AU


    ############################################################
    # Quick check?
    fig=plt.figure()
    ax=plt.axes()
    #ax.plot(r_au,brightness/max(brightness),'*',markersize=0.9,label="model")
    #if mode=='jband':
    #    ax.plot(r_au,brightness,'*',markersize=0.9,label="model")
    ax.set_xlabel(r"$r$ (AU)")
    if mode=='alma':
        ax.plot(r_au,brightness/max(brightness),'*',markersize=0.9,label="model")
        odata=np.loadtxt("alma_profile_observed_normed.dat")
        r_obs=odata[:,0:1]
        b_obs=odata[:,1:2]
        ax.plot(r_obs,b_obs,label="observation")
        ax.set_title("No convolved")
        ax.axvline(74)
        ax.legend()
    elif mode=='jband':
        factor=1/50.8849
        brightness=brightness*factor
        ax.plot(r_au,brightness,'*',markersize=0.9,label="model")
        odata=np.loadtxt("jband_profile_observed_normed.dat")
        r_obs=odata[:,0:1]
        b_obs=odata[:,1:2]
        ax.axvline(54)
        ax.plot(r_obs,b_obs,label='observation')
        ax.set_ylim(-0.1,1.8)
        ax.legend()
    plt.show()

    return None

alma_data('/data/users/bportilla/runs/runs/run_21/','RTout0001_000855.84.fits.gz','alma',150)
alma_data('/data/users/bportilla/runs/runs/run_21/','RTout0001_000001.25.fits.gz','jband',150)
