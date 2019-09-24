import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl
plt.style.use('fancy')

############################################################
# Canonical units
uL=6.9551e8
uM=1.989e30
uT=2.639e35

############################################################
# Constants and definitions 
AU=149.6e9
pc=206265.0*AU
ly=63241.1*AU
c=300.0e6 
k=1.3806e-23 
h=6.6261e-34
R=1.26*6.9551e8
dist=113.43*pc
Teff=3972.0

############################################################
# MCMax3D data
data=np.loadtxt('star0001.dat')

flux_min=[]
for i in range(0,data.shape[0]):
    factor=c/((data[i][0]*1e-6)**2) # (m/s/m^2)
    flux_min.append(factor*(data[i][1]*1e-26)) # from Jy to W/m^3
flux_min=np.array(flux_min)

############################################################
# Wavelenght range (input in microns)
lmin=0.1
lmax=1e3
Npoints=1e5
lrange=np.linspace(lmin/1e6,lmax/1e6,Npoints) # In mts

############################################################
# Convertion to canonical units
c=c*uT/uL
k=k*uT**2/(uL**2*uM)
h=h*uT/(uL**2*uM)
R=R/uL
dist=dist/uL
lrange=lrange/uL

############################################################
# Planck function
def planck_function(x,Teff): 
    value=2*np.pi*h*c**2/x**5 * 1/(np.exp(h*c/(x*k*Teff))-1) 
    return value

############################################################
# Flux at r=dist
flux=planck_function(lrange,Teff)*(R/dist)**2

############################################################
# Converting to SI
flux=flux*uM/(uL*uT**3) # W/m^3
lrange=lrange*uL # m

############################################################
# Converting wavelenght to micron and flux to W/m^2/micron
flux=flux*1e-6 # W/m^2/micron
lrange=lrange*1e6 # microns

############################################################
# Generating file of theoretical spectrum for MCMax3D
file=open('black_body_PDS70.dat','w')
for i in range(0,len(lrange)):
    #file.write("%.16e %.16e\n"%(lrange[i],(flux*lrange)[i]))
    file.write("%.16e %.16e\n"%(lrange[i],flux[i]))


############################################################
# Plot lambda*F_lambda vs wavelenght

fig,(ax_1,ax_2) = plt.subplots(1, 2,figsize=(16,8),sharex=True)

ax_1.plot(lrange,flux*lrange)
ax_1.set_xscale('log')
ax_1.set_yscale('log')
ax_1.set_xlabel(r'$\lambda \, (\mu\mathrm{m})$')
ax_1.set_ylabel(r'$\lambda F_{\lambda}$ (W/m^2)')
ax_1.set_ylim(1e-17,1e-12)

ax_2.plot(lrange,flux,label='my spectrum')
ax_2.plot(data[:,0:1],flux_min/1e6,'.',label='MCMax3D')

ax_2.legend()
ax_2.set_xscale('log')
ax_2.set_yscale('log')
ax_2.set_xlabel(r'$\lambda \, (\mu\mathrm{m})$')
ax_2.set_ylabel(r'$F_{\lambda}$ (W/m^2/\mu\mathrm{m})')
#ax_2.set_ylim(1e-17,1e-12)

plt.show()
