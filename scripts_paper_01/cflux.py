import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits
from photutils import EllipticalAnnulus,CircularAnnulus,EllipticalAperture
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_image
import sys
plt.style.use('fancy')

############################################################
# Load MCMax3D image
file="/data/users/bportilla/runs/runs/run_23_copy/output/RTout0001_000855.84.fits.gz"
pxsize=4/1000.0 # arcsec/px
theta=49.7
hdulist=fits.open(file)
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
# Get radial profiles

## observation
file="../observations/PDS70_cont-final.fits"
pxsize=20.0
pa=158.6
theta=49.7
hdulist=fits.open(file)
data_obs=hdulist[0].data[0][0]*1000 # mJy

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
a_in_array=[]
for i in np.arange(yc+dr,yc+0.5*pixel_lim,dr):
    a_in_array.append(i-xc)
a_out_array=[i+dr for i in a_in_array]
b_out_array=[i*(1-e**2)**0.5 for i in a_out_array]

apertures=[EllipticalAnnulus((yc,xc),a_in=ain,a_out=aout,b_out=bout,theta=angle_annulus)
           for (ain,aout,bout) in zip(a_in_array,a_out_array,b_out_array)]


############################################################
# Do a check
aperture=apertures[-5]
plt.imshow(data_obs)
plt.title("Image observation")
aperture.plot(color='red',lw=1)
plt.show()

# Radial distance of each annulus
r_arcsec=[(j+0.5*(i-j))*pxsize for (i,j) in zip(a_out_array,a_in_array)] # arcsec
r_rad=[(i*units.arcsec).to(units.rad).value for i in r_arcsec] # rad
r_au=[((i*d)*units.pc).to(units.au).value for i in r_rad] # AU

# Creating numpy arrays
r_au=np.array(r_au)
r_arcsec=np.array(r_arcsec)

phot_table=aperture_photometry(data_obs,apertures)
col_values=[]
for col in phot_table.colnames:
    col_values.append(phot_table[col][0])
brightness=[col_values[i] for i in range(3,len(col_values))]


############################################################
# Creating brightness profile
fig=plt.figure()
ax=plt.axes()
ax.plot(r_au,brightness/max(brightness),'*')
ax.set_xlabel(r"Projected radial distance (AU)")
ax.set_ylabel("Density flux (mJy/beam)")
ax.set_title("Radial profile observation")
plt.show()

############################################################
# Creating file
file=open('alma_radial_profile_observed.dat',"w")
for i in range(0,len(r_au)):
    file.write('%.15e %.15e \n'%(r_au[i],brightness[i]/max(brightness)))


############################################################
############################################################
############################################################
############################################################

## model
pxsize=4.0/1000.0 # arcsec/pixel
pa=158.6
theta=49.7

############################################################
# Input params
pa=158.6
angle_annulus=0.0#((pa-90.0)*units.deg).to(units.rad).value # Position angle
#angle_annulus=((pa)*units.deg).to(units.rad).value # Position angle
e=np.sin((theta*units.deg).to(units.rad).value) # eccentricity of the annulus

# Determining limit for radial profile
linear_lim=2*(120) # AU
angular_lim=linear_lim/((d*units.pc).to(units.au).value) # rad
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
aperture=apertures[-5]
plt.imshow(data_mod)
plt.title("Image model")
aperture.plot(color='red',lw=1)
plt.show()

# Radial distance of each annulus
r_arcsec=[(j+0.5*(i-j))*pxsize for (i,j) in zip(a_out_array,a_in_array)] # arcsec
r_rad=[(i*units.arcsec).to(units.rad).value for i in r_arcsec] # rad
r_au=[((i*d)*units.pc).to(units.au).value for i in r_rad] # AU

# Creating numpy arrays
r_au=np.array(r_au)
r_arcsec=np.array(r_arcsec)

phot_table=aperture_photometry(data_mod,apertures)
col_values=[]
for col in phot_table.colnames:
    col_values.append(phot_table[col][0])
brightness=[col_values[i] for i in range(3,len(col_values))]


############################################################
# Creating brightness profile
fig=plt.figure()
ax=plt.axes()
ax.plot(r_au,brightness/max(brightness),'*')
ax.set_xlabel(r"Projected radial distance (AU)")
ax.set_ylabel("Density flux (mJy/beam)")
ax.set_title("Radial profile model")
plt.show()


############################################################
# Creating file
file=open('alma_radial_profile_modeled.dat',"w")
for i in range(0,len(r_au)):
    file.write('%.15e %.15e \n'%(r_au[i],brightness[i]/max(brightness)))




