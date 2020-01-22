import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits
from photutils import EllipticalAnnulus,CircularAnnulus,EllipticalAperture
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_image
import sys
plt.style.use('fancy')

def alma():
    lim=120.0

    ############################################################
    # Absolute paths to files
    fits_image='RTout0001_000855.84.fits.gz'
    path_fits_image='output/'+fits_image
    path_image_file='Image_alma.out'
    path_input_file='input.dat'

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
    data_mod=hdulist[0].data[0]

    ############################################################
    # mJy/arcsec2 -> mJy/beam
    beamx=74.0/1000 # arcsec
    beamy=57.0/1000 # arcsec
    beam_size=np.pi*beamx*beamy # arcsec2
    data_mod=data_mod*beam_size

    ############################################################
    # Convolve?
    data_mod=convolve_image(data_mod)

    ############################################################
    #
    # Get radial profiles - model
    #
    ############################################################

    ############################################################
    # Input params
    angle_annulus=0.0#((pa-90.0)*units.deg).to(units.rad).value # Position angle

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

    """
    ############################################################
    # Do a check
    a=0.01
    vmin_jband=np.percentile(data_mod,a)
    vmax_jband=np.percentile(data_mod,100-a)
    aperture=apertures[-1]
    plt.imshow(data_mod,clim=(vmin,vmax))
    plt.title("Image model")
    aperture.plot(color='red',lw=1)
    plt.show()
    """

    phot_table=aperture_photometry(data_mod,apertures)
    col_values=[]
    for col in phot_table.colnames:
        col_values.append(phot_table[col][0])
    brightness=[col_values[i] for i in range(3,len(col_values))]

    # Radial distance of each annulus
    r_arcsec=[(j+0.5*(i-j))*pxsize for (i,j) in zip(a_out_array,a_in_array)] # arcsec
    r_rad=[(i*units.arcsec).to(units.rad).value for i in r_arcsec] # rad
    r_au=[(i*d) for i in r_rad] # AU

    # Creating numpy arrays
    r_au=np.array(r_au)
    r_arcsec=np.array(r_arcsec)


    """
    ############################################################
    # Creating brightness profile
    fig=plt.figure()
    ax=plt.axes()
    ax.plot(r_au,brightness/max(brightness),'*')
    ax.set_xlabel(r"Projected radial distance (AU)")
    ax.set_ylabel("Density flux (mJy/beam)")
    ax.set_title("Radial profile model")
    plt.show()
    """

    ############################################################
    # Creating file
    file=open('alma_radial_profile_modeled.dat',"w")
    for i in range(0,len(r_au)):
        file.write('%.15e %.15e \n'%(r_au[i],brightness[i]/max(brightness)))


    return data_mod

alma()
