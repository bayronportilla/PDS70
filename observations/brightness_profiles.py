import numpy as np
import sys
from astropy import units
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
import matplotlib.ticker as ticker
from photutils import EllipticalAnnulus
from photutils import aperture_photometry
#plt.style.use('fancy')


############################################################
# Fancy style
plt.rc('font', **{'family': 'family', 'family': ['serif']})
plt.rc('text', usetex=True) #uses Latex instead of Tex to compile axes labels

def radial_profile_observation(file,mode):

    # Loading FITS image
    hdulist=fits.open(file)
    #data=(hdulist[0].data)[0][0]
    data=(hdulist[0].data)

    #plt.imshow(data,origin="lower")
    #plt.show()
    #sys.exit()

    """
    # Removing NaN
    ave_flux=1e-5
    for i in range(0,data.shape[0]):
        for j in range(0,data.shape[1]):
            if np.isnan(data[i][j])==True:
                data[i][j]=ave_flux
            else:
                continue
    """
    
    ############################################################
    # Input params
    #pixel_scale=0.02 # arcsec per pixel
    pixel_scale=0.01226 # arcsec per pixel
    d=113.43 # Distance 
    PA=158.6 # Position angle
    #e=0.775 # eccentricity of the annulus
    e=0.767 # eccentricity of the annulus
    
    ############################################################
    # Creating image
    x_array=np.linspace(0,data.shape[0],data.shape[0])
    y_array=np.linspace(0,data.shape[1],data.shape[1])
    x_array=[i*pixel_scale for i in x_array]
    y_array=[i*pixel_scale for i in y_array]
    gauge_x=0.5*max(x_array)
    gauge_y=0.5*max(y_array)
    x_array=[i-gauge_x for i in x_array]
    y_array=[i-gauge_y for i in y_array]
    extent_angular=[max(x_array),min(x_array),min(y_array),max(y_array)]
    fig=plt.figure()
    ax=plt.axes()
    image_1=ax.imshow(1000*data,cmap='hot',origin='lower',extent=extent_angular)
    cbar=fig.colorbar(image_1,pad=0.0,aspect=25)
    #cbar.ax.set_ylabel('Intensity (mJy/beam)')
    #ax.set_xlim(1.5,-1.5)
    #ax.set_ylim(-1.5,1.5)
    cbar.ax.set_ylabel('Intensity (???)')
    ax.set_xlim(1,-1)
    ax.set_ylim(-1,1)
    ax.set_xlabel(r'$\Delta \mathrm{RA}$ (arcsec)')
    ax.set_ylabel(r'$\Delta \mathrm{Dec}$ (arcsec)')
    fig.savefig('%s_image.png'%file)

    
    ############################################################
    # Creating brightness profile
    theta_phot=((PA+90)*units.deg).to(units.rad) # angle for photometry measured form x-axis
    x_center=0.5*data.shape[0] # Image center in data coordinates
    y_center=0.5*data.shape[1] # Image center in data coordinates
    dr=1 # Width of the annulus
    ain_array=[] # array of inner semi-major axis
    r_array=[] # array of (linear) radial distance
    for i in np.arange(x_center+dr,data.shape[0],dr): 
        ain_array.append(i-x_center)
    aout_array=[i+dr for i in ain_array]
    bout_array=[i*(1-e**2)**0.5 for i in aout_array]
   
    # Creating apertures
    apertures=[EllipticalAnnulus((x_center,y_center),a_in=ain,a_out=aout,b_out=bout,theta=theta_phot.value)
               for (ain,aout,bout) in zip(ain_array,aout_array,bout_array)]

    r_array=[((i-j)*0.5+j)*pixel_scale for i,j in zip(aout_array,ain_array)]
    r_array=[(i*units.arcsec).to(units.rad).value for i in r_array]
    r_array=[i*((d*units.pc).to(units.au).value) for i in r_array]

    phot_table=aperture_photometry(data,apertures)
    col_values=[]
    for col in phot_table.colnames:
        col_values.append(phot_table[col][0])
    brightness=[col_values[i] for i in range(3,len(col_values))]

   
    ############################################################
    # Creating brightness profile
    fig=plt.figure()
    ax=plt.axes()
    ax.plot(r_array,brightness)
    ax.set_xlabel(r"$r$ (AU)")
    ax.set_ylabel("Intensity (???)")
    #ax.set_xlim(0,150)
    ax.set_xlim(0,110)
    fig.savefig('%s_radial_profile.png'%file)
          
    return "Files generated!"

#radial_profile_observation("PDS70_cont-final.fits",'Q')
radial_profile_observation("PDS_70_2016-03-26_QPHI_CROPPED.fits",'Q')
    

