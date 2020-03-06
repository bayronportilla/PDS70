import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits
from photutils import EllipticalAnnulus,CircularAnnulus,EllipticalAperture,RectangularAperture
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_image
import matplotlib.gridspec as gridspec
import sys
from astropy.table import Table
plt.style.use('fancy')

def toAU(x,y,xc,yc,pxsize,d):
    dr=((x-xc)**2+(y-yc)**2)**0.5
    dr=(dr*pxsize)
    dr=(dr*units.arcsec).to(units.rad).value
    dr=dr*((d*units.pc).to(units.au).value)
    return dr

def toPX(l,pxsize,d):
    dl=(l*units.au).to(units.pc).value
    dl=((dl/d)*units.rad).to(units.arcsec).value
    dl=dl/pxsize
    return dl

############################################################
# Get radial profiles of the observation
file="../observations/PDS_70_2017-08-01_QPHI_amorph.fits"
pxsize=12.26
pa=158.6
theta=49.7
hdulist=fits.open(file)
data_obs=hdulist[0].data # adu's

############################################################
# Input params
pxsize=pxsize/1000.0 # arcsecs per pixel
d=113.43 # Distance in pc 
pa=158.6
angle_annulus=((pa-90.0)*units.deg).to(units.rad).value # Position angle
#angle_annulus=((pa)*units.deg).to(units.rad).value # Position angle
e=np.sin((theta*units.deg).to(units.rad).value) # eccentricity of the annulus

# Determining limit for radial profile
linear_lim=2*(120) # AU
angular_lim=linear_lim/((d*units.pc).to(units.au).value) # rad
angular_lim=(angular_lim*units.rad).to(units.arcsec).value # arcsec
pixel_lim=int(round(angular_lim/pxsize))

xc=0.5*data_obs.shape[0] # Image center in data coordinates
yc=0.5*data_obs.shape[1] # Image center in data coordinates
dr=1.0 # Width of the annulus
w=1.0
h=1.0
lim=100
lim=toPX(lim,pxsize,d)
InRad=np.pi/180.0
xc_array=[]
yc_array=[]

x0=xc
y0=yc
xc_array.append(x0)
yc_array.append(y0)

for l in np.arange(1,lim,1):
    
    xval=x0+l*np.cos(angle_annulus+1.0*np.pi)
    yval=y0+l*np.sin(angle_annulus+1.0*np.pi)
    xc_array.append(xval)
    yc_array.append(yval)
    
positions=[(i,j) for (i,j) in zip(xc_array,yc_array)]

InRad=np.pi/180.0
#apertures=RectangularAperture(positions,w,h,pa*InRad) 
apertures=RectangularAperture(positions,w,h,angle_annulus) 
#apertures=[RectangularAperture((i,j),w,h,pa*InRad) for (i,j) in zip(xc_array+xc,yc_array+yc) ]
#aperture=RectangularAperture((yc,xc),w,h,angle_annulus)

############################################################
# Do a check
a=0.01
vmin=np.percentile(data_obs,a)
vmax=np.percentile(data_obs,100-a)
plt.imshow(data_obs,clim=(vmin,vmax))
plt.xlim()
plt.title("Image observation")
apertures.plot(color='red',lw=1)
plt.show()


phot_table=aperture_photometry(data_obs,apertures)
r_au=[toAU(phot_table['xcenter'][i].value,phot_table['ycenter'][i].value,xc,yc,pxsize,d) for i in range(0,len(phot_table))]
brightness=[phot_table['aperture_sum'][i] for i in range(0,len(phot_table))]
brightness=np.array(brightness)


############################################################
# Creating brightness profile normalized
fac=1/3.41577
brightness=brightness*fac
fig=plt.figure()
ax=plt.axes()
ax.plot(r_au,brightness,'-')
ax.set_xlabel(r"Projected radial distance (AU)")
ax.set_ylabel("Density flux (a.u.)")
ax.set_title("Radial profile observation")
ax.set_ylim(-0.1,5)
plt.show()

sys.exit()
############################################################
# Creating file
file=open('jband_radial_profile_observed_HR.dat',"w")
for i in range(0,len(r_au)):
    file.write('%.15e %.15e \n'%(r_au[i],brightness[i]))
