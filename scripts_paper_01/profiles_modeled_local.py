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
from cflux_alma_local import alma
from cflux_jband_local import jband
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

   
    data_alma=alma()
    data_jband=jband()
    
    ############################################################
    # Dynamic range
    a_alma=0.01
    vmin_alma=np.percentile(data_alma,a_alma)
    vmax_alma=np.percentile(data_alma,100-a_alma)

    a_jband=0.01
    vmin_jband=np.percentile(data_jband,a_jband)
    vmax_jband=np.percentile(data_jband,100-a_jband)
    
  
    ############################################################
    # Observed profile

    odata=np.loadtxt("alma_radial_profile_observed.dat")
    r_obs=odata[:,0:1]
    b_obs=odata[:,1:2]

    odata_j=np.loadtxt("jband_radial_profile_observed.dat")
    r_obs_j=odata_j[:,0:1]
    b_obs_j=odata_j[:,1:2]

    ############################################################
    # Loading profiles
    mprofile_alma=np.loadtxt("alma_radial_profile_modeled.dat")
    r_alma=np.reshape(mprofile_alma[:,0:1],mprofile_alma.shape[0])
    b_alma=np.reshape(mprofile_alma[:,1:2],mprofile_alma.shape[0])
    
    mprofile_jband=np.loadtxt("jband_radial_profile_modeled.dat")
    r_jband=np.reshape(mprofile_jband[:,0:1],mprofile_jband.shape[0])
    b_jband=np.reshape(mprofile_jband[:,1:2],mprofile_jband.shape[0])

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
    ax3.plot(r_alma,b_alma,'*',markersize=0.9,label="model")
    ax4.plot(r_jband,b_jband,'*',markersize=0.9,label="model")
    ax4.plot(r_obs_j,b_obs_j,'*',markersize=0.9,label="observation")
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
    ax4.set_ylim(-0.11,2.0)

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

print(output_data())

