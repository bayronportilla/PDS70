import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits
from photutils import EllipticalAnnulus,CircularAnnulus,EllipticalAperture
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_model
import sys
plt.style.use('fancy')

def topx(l,pxsize,d):
    x=(l*units.au).to(units.pc).value
    x=((x/d)*units.rad).to(units.arcsec).value
    x=x/pxsize
    return x


def get_profile(image,pxsize,amean,width,inc,PA,d):

    ############################################################
    #
    # Inputs.
    # image: fits file of the image.
    # pxsize: pixel scale of the image in arcsex/px.
    # amean: mean semi-major axis of the elliptical aperture (AU).
    # width: total width of the aperture (AU).
    # inc: inclination of the disk (deg).
    # PA: east-north measured position angle of the disk (deg).
    # d: distance to the source in pc
    #
    # Returns 
    # 
    # 
    #
    ############################################################


    ############################################################
    # Loading data
    hdu=fits.open(image)
    data=hdu[0].data[0][0]

    
    ############################################################
    # Derived quantities
    x0=data.shape[0]*0.5
    y0=data.shape[1]*0.5
    inc=(inc*units.deg).to(units.rad).value
    e=np.sin(inc)


    ############################################################
    # Creating elliptical aperture
    amin=amean-width*0.5 # AU
    amax=amean+width*0.5 # AU
    bmax=amax*(1-e**2)**0.5 # AU

    amin=topx(amin,pxsize,d) # px
    amax=topx(amax,pxsize,d) # px
    bmax=topx(bmax,pxsize,d) # px

    angle=((PA+90)*units.deg).to(units.rad).value

    aperture=EllipticalAnnulus((x0,y0),a_in=amin,a_out=amax,b_out=bmax,theta=angle)
    
    # Do a check?
    plt.imshow(data)
    aperture.plot(color='red',lw=1)
    plt.show()
    
    

    ############################################################
    # Creating aperture mask
    mask=aperture.to_mask(method="center")
    
    # Do a check?
    plt.imshow(mask)
    plt.colorbar()
    plt.show()
    

    ############################################################
    # Extracting pixels located inside the aperture
    aperture_data=mask.multiply(data)
    print(aperture_data.shape)
    # Do a check?
    plt.imshow(aperture_data)
    plt.colorbar()
    plt.show()
    

#    print(mask)


    
    #plt.imshow(data)
    #plt.show()

    return 0

get_profile("../observations/PDS70_cont-final.fits",0.020,74.0,10.0,49.7,158.6,113.43)
