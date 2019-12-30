import numpy as np
import sys
from astropy import units
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
import matplotlib.ticker as ticker
from photutils import EllipticalAnnulus
from photutils import aperture_photometry
plt.style.use('fancy')

def radial_profile_observation(file,mode,theta,PA):

    # Loading FITS image
    hdulist=fits.open(file)
    data=(hdulist[0].data)[0][0]

    # Removing NaN
    ave_flux=1e-5
    for i in range(0,data.shape[0]):
        for j in range(0,data.shape[1]):
            if np.isnan(data[i][j])==True:
                data[i][j]=ave_flux
            else:
                continue
    
    ############################################################
    # Input params
    pxsize=20.0/1000 # arcsecs per pixel
    d=113.43 # Distance in pc 
    angle_annulus=((PA-90.0)*units.deg).to(units.rad).value # Position angle
    #angle_annulus=((PA)*units.deg).to(units.rad).value # Position angle
    e=np.sin((theta*units.deg).to(units.rad).value) # eccentricity of the annulus
    
    ############################################################
    # Creating image
    """
    x_array=np.linspace(0,data.shape[0],data.shape[0])
    y_array=np.linspace(0,data.shape[1],data.shape[1])
    x_array=[i*pxsize for i in x_array]
    y_array=[i*pxsize for i in y_array]
    gauge_x=0.5*max(x_array)
    gauge_y=0.5*max(y_array)
    x_array=[i-gauge_x for i in x_array]
    y_array=[i-gauge_y for i in y_array]
    """

    fig=plt.figure()
    ax=plt.axes()
    fov=data.shape[0]*pxsize # arcsec
    extent_angular=[0.5*fov,-0.5*fov,0-0.5*fov,0.5*fov]
    a=0.01
    vmin=np.percentile(data,a)
    vmax=np.percentile(data,100-a)
    image_1=ax.imshow(data,origin='lower',cmap='hot',extent=extent_angular,clim=(vmin,vmax))
    cbar=fig.colorbar(image_1,pad=0.0,aspect=25)
    cbar.ax.set_ylabel('Intensity (Jy/beam)')

    # Zoom
    lim=1.5
    ax.set_xlim(lim,-lim)
    ax.set_ylim(-lim,lim)

    ax.set_xlabel(r'$\Delta \mathrm{RA}$ (arcsec)')
    ax.set_ylabel(r'$\Delta \mathrm{Dec}$ (arcsec)')
    fig.savefig('image_observed.png')
    
    
    ############################################################
    # Creating brightness profile

    # Determining limit for radial profile
    linear_lim=2*(150) # AU
    angular_lim=linear_lim/((d*units.pc).to(units.au).value) # rad
    angular_lim=(angular_lim*units.rad).to(units.arcsec).value # arcsec
    pixel_lim=int(angular_lim/pxsize)
    
    xc=0.5*data.shape[0] # Image center in data coordinates
    yc=0.5*data.shape[1] # Image center in data coordinates
    dr=0.1 # Width of the annulus
    ain_array=[] # array of inner semi-major axis
    r_array=[] # array of (linear) radial distance
    
    #return (xc,pixel_lim,data.shape[0])
    #sys.exit()
    for i in np.arange(xc+dr,xc+0.5*pixel_lim,dr): 
    #for i in np.arange(xc+dr,data.shape[0],dr): 
        ain_array.append(i-xc)
    aout_array=[i+dr for i in ain_array]
    bout_array=[i*(1-e**2)**0.5 for i in aout_array]
   
    # Radial distance of each annulus
    r_arcsec=[(j+0.5*(i-j))*pxsize for (i,j) in zip(aout_array,ain_array)] # arcsec
    r_rad=[(i*units.arcsec).to(units.rad).value for i in r_arcsec] # rad
    r_au=[((i*d)*units.pc).to(units.au).value for i in r_rad] # AU

    # Creating numpy arrays
    ain_array=np.array(ain_array)
    aout_array=np.array(aout_array)
    bout_array=np.array(bout_array)
    r_au=np.array(r_au)
    r_arcsec=np.array(r_arcsec)

    # Creating elliptical annulus apertures
    apertures=[EllipticalAnnulus((xc,yc),a_in=ain,a_out=aout,b_out=bout,theta=angle_annulus)
               for (ain,aout,bout) in zip(ain_array,aout_array,bout_array)]

    phot_table=aperture_photometry(data,apertures)
    col_values=[]
    for col in phot_table.colnames:
        col_values.append(phot_table[col][0])
    brightness=[col_values[i] for i in range(3,len(col_values))]

   
    ############################################################
    # Creating brightness profile
    fig=plt.figure()
    ax=plt.axes()
    ax.plot(r_au,brightness,'-')
    ax.set_xlabel(r"$r$ (AU)")
    ax.set_ylabel("Intensity")
    
    #ax.set_ylim(0,3)

    fig.savefig('radial_profile.png')

    norm_factor=1/max(brightness)
    brightness_normed=[i*norm_factor for i in brightness]

    file=open('radial_profile_observed.dat',"w")
    for i in range(0,len(r_au)):
        #file.write('%.15e %.15e \n'%(r_au[i],brightness_normed[i]))
        file.write('%.15e %.15e \n'%(r_au[i],brightness[i]))

    return "Files generated!"

print(radial_profile_observation("PDS70_cont-final.fits",'Q',49.7,158.6))
#radial_profile_observation("PDS_70_2016-03-26_QPHI_CROPPED.fits",'Q')
    
