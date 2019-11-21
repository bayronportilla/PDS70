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

def radial_profile_observation(file,mode):

    # Loading FITS image
    hdulist=fits.open(file)    
    data=(hdulist[0].data)[0][0]
        
    #a=0.15
    #vmin=np.percentile(data,a)
    #vmax=np.percentile(data,100-a)
    
    # Removing NaN
    ave_flux=1e-5
    for i in range(0,data.shape[0]):
        for j in range(0,data.shape[1]):
            if np.isnan(data[i][j])==True:
                data[i][j]=ave_flux
            else:
                continue

    # Brightness profile
    pixel_scale=0.01226 # arcsec per pixel
    d=113.43 # Distance 
    PA=158.6 # Position angle
    theta_phot=((PA+90)*units.deg).to(units.rad) # angle for photometry measured form x-axis
    e=0.775 # eccentricity of the annulus
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
    sys.exit()
    # Creating figures
    #fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,5))
    fig=plt.figure()
    ax=plt.axes()
    image_1=ax.imshow(1000*data,cmap='hot',origin='lower',extent=extent_angular)
    plt.show()
    sys.exit()
    ax1.set_xlabel(r'$\Delta \mathrm{RA}$ (arcsec)')
    ax1.set_ylabel(r'$\Delta \mathrm{Dec}$ (arcsec)')
            
    # Color bar
    #cbar=fig.colorbar(image_1,pad=0.0,aspect=25)
    #cbar.ax1.set_ylabel('Intensity (mJy/beam)')

    #plt.xlim(1.5,-1.5)
    #plt.ylim(-1.5,1.5)
    #plt.show()
    #sys.exit()
    #ax.plot(r_au,brightness,'.',markersize=0.8)
    #ax.set_xlabel(r"$r$ (AU)")
    #ax.set_ylabel("Intensity (???)")
    
    #ax2=ax.twiny()
    #ax2.plot(r_arcsec,brightness)
    #ax2.set_xlabel(r"$r$ (arcsec)")
    fig.savefig('prueba.png')
    
    #plt.show()
    
    sys.exit()
    """
    for i in phot_table.colnames:
        print(phot_table['aperture_sum_%d'%index][0])
        index+=1

    #print(phot_table)

    # Doing aperture photometry
    #phot_table=aperture_photometry(data,aperture)
    #phot_table['aperture_sum'].info.format = '%.8g'
    #print(phot_table)
    

    #plt.show()


    #fig.savefig('RTout0001_000863.31_%s.png'%(mode))
    #fig.savefig('RTout0001_000863.31_%s_%d.png'%(mode,int(phi)))
    #fig.savefig('lower.png')
          
    return "Image generated!"
    """
radial_profile_observation("PDS70_cont-final.fits",'Q')
#radial_profile_observation("PDS_70_2016-03-26_QPHI_CROPPED.fits",'Q')
    
