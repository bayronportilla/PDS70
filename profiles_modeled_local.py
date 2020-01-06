import numpy as np
import sys
from astropy import units
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
from photutils import EllipticalAnnulus
from mcmax3d_analysis.mcmax3d_observables import convert_flux
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_image
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
plt.style.use('fancy')

def output_data():
    
    ############################################################
    # Converting MCMax3D spectra 
    convert_flux('output/MCSpec0001.dat','spectrum_PDS70_system.dat')
    convert_flux('output/star0001.dat','spectrum_PDS70_star.dat')


    ############################################################
    # Loading surface density profile
    surface_density=np.loadtxt('surface_density_PDS70.dat')


    ############################################################
    # Loading data (converted) from MCMax3D 
    system_spectrum=np.loadtxt('spectrum_PDS70_system.dat')
    stellar_spectrum=np.loadtxt('spectrum_PDS70_star.dat')

    ############################################################
    # Loading photometric data
    photometric_data=np.loadtxt('PDS70_photometry.dat')
    x_photometric=np.reshape(photometric_data[:,0:1],photometric_data.shape[0])
    y_photometric=np.reshape(photometric_data[:,1:2],photometric_data.shape[0])
    y_photometric_error=np.reshape(photometric_data[:,2:3],photometric_data.shape[0])

    ############################################################
    # Creating arrays for plotting
    x_system=np.reshape(system_spectrum[:,0:1],system_spectrum.shape[0])
    y_system=np.reshape(system_spectrum[:,1:2],system_spectrum.shape[0])
    x_star=np.reshape(stellar_spectrum[:,0:1],stellar_spectrum.shape[0])
    y_star=np.reshape(stellar_spectrum[:,1:2],stellar_spectrum.shape[0])

    ############################################################
    # Creating array of errors
    y_system_error=x_photometric*y_photometric_error

    lim=120.0

    ############################################################
    # Absolute paths to files
    fits_alma='RTout0001_000855.84.fits.gz'
    fits_jband='RTout0001_000001.25.fits.gz'
    path_fits_file_alma='output/'+fits_alma
    path_fits_file_jband='output/'+fits_jband
    path_image_file='Image.out'
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

    ############################################################
    # For photometry
    e=np.sin(theta) # eccentricity of the annulus

    ############################################################
    # Load data
    hdulist_alma=fits.open(path_fits_file_alma)
    data_alma=hdulist_alma[0].data[0]
    data_alma=convolve_image(data_alma)
    
    hdulist_jband=fits.open(path_fits_file_jband)
    data_Q=hdulist_jband[0].data[1]
    data_U=hdulist_jband[0].data[2]
    data_jband=data_Q*np.cos(2*phi) + data_U*np.sin(2*phi)

    ############################################################
    # Brightness profile
    linear_lim=2*(lim) # AU
    angular_lim=((linear_lim/d)*units.rad).to(units.arcsec).value # arcsec
    pixel_lim=int(round(angular_lim/pxsize))

    # Creating elliptical aperture object in pixel coordinates
    xc=0.5*data_alma.shape[0]
    yc=0.5*data_alma.shape[1]
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
    a_alma=0.01
    vmin_alma=np.percentile(data_alma,a_alma)
    vmax_alma=np.percentile(data_alma,100-a_alma)

    a_jband=0.01
    vmin_jband=np.percentile(data_jband,a_jband)
    vmax_jband=np.percentile(data_jband,100-a_jband)

    """
    ############################################################
    # Do a check?
    aperture=apertures[-1]
    fig=plt.figure()
    ax=plt.axes()
    ax.imshow(data_alma,clim=(vmin,vmax))
    aperture.plot(color='red',lw=1)
    """

    ############################################################
    # Photometry for alma
    phot_table_alma=aperture_photometry(data_alma,apertures)
    col_values_alma=[]
    for col in phot_table_alma.colnames:
        col_values_alma.append(phot_table_alma[col][0])
    brightness_alma=[col_values_alma[i] for i in range(3,len(col_values_alma))]

    ############################################################
    # Photometry for jband
    phot_table_jband=aperture_photometry(data_jband,apertures)
    col_values_jband=[]
    for col in phot_table_jband.colnames:
        col_values_jband.append(phot_table_jband[col][0])
    brightness_jband=[col_values_jband[i] for i in range(3,len(col_values_jband))]

    # Create arrays for output file
    r_arcsec=[(j+0.5*(i-j))*pxsize for (i,j) in zip(aout_array,ain_array)] # arcsec
    r_rad=[(i*units.arcsec).to(units.rad).value for i in r_arcsec] # rad
    r_au=[((i*d)*units.au).value for i in r_rad] # AU

    ############################################################
    # Observed profile

    odata=np.loadtxt("alma_radial_profile_observed.dat")
    r_obs=odata[:,0:1]
    b_obs=odata[:,1:2]

    ############################################################
    # Plotting
    fig=plt.figure()
    gs=gridspec.GridSpec(2,3)
    ax1=plt.subplot(gs[0,0:1])
    ax2=plt.subplot(gs[0,1:2])
    ax3=plt.subplot(gs[0,2:3])
    ax4=plt.subplot(gs[1,0:1])
    ax5=plt.subplot(gs[1,1:2])
    ax6=plt.subplot(gs[1,2:3])
    ax1.imshow(data_alma,clim=(vmin_alma,vmax_alma))
    ax2.imshow(data_jband,clim=(vmin_jband,vmax_jband))
    ax3.plot(r_obs,b_obs,label="observation")
    ax3.plot(r_au,brightness_alma/max(brightness_alma),'*',markersize=0.9,label="model")
    ax4.plot(r_au,brightness_jband/max(brightness_jband),'*',markersize=0.9,label="model")
    ax4.axvline(54.0)
    ax5.plot(x_system,x_system*y_system,label='MCMax3D output',linewidth=1.5)
    ax5.plot(x_star,x_star*y_star,color='orange')
    ax5.errorbar(x_photometric,x_photometric*y_photometric,yerr=y_system_error,label='Photometric data',fmt='.',markersize=7,color='green')
    ax6.plot(surface_density[:,0:1],surface_density[:,1:2])

    ############################################################
    # Scales
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax6.set_xscale('log')
    ax6.set_yscale('log')

    ############################################################
    # lims
    ax5.set_ylim(1e-17,1e-12)
    ax6.set_ylim(1e-6,1e5)

    ############################################################
    # Legend
    ax3.legend()
    ax4.legend()
    
    ############################################################
    # Labels
    ax3.set_xlabel(r"$r$ (AU)")
    ax3.set_ylabel(r"Normalized density flux")
    ax4.set_xlabel(r"$r$ (AU)")
    ax4.set_ylabel(r"Normalized density flux")
    ax5.set_xlabel(r'$\lambda \, (\mu\mathrm{m})$')
    ax5.set_ylabel(r'$\lambda F_{\lambda}$ (W/m^2)')
    ax6.set_xlabel(r'$r$ (AU)')
    ax6.set_ylabel(r'$\Sigma_{\mathrm{dust}}$ (g/cm^2)')
    plt.show()
    
    fig.savefig("fig_all.png")

    return None

output_data()

