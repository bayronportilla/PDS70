import numpy as np
import matplotlib.pyplot as plt
import sys

############################################################
# Canonical units
uL=6.9551e8
uM=1.989e30
uT=2.639e35

############################################################
# Constants and definitions 
c=300.0e6 
k=1.3806e-23 
h=6.6261e-34
R=1.26*6.9551e8
dist=23396638.95*149.6e9
Teff=3972.0

############################################################
# Wavelenght range (microns)
lmin=0.1
lmax=1e3
Npoints=1e6
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
# Flux in W/m^3
flux=planck_function(lrange,Teff)*(R/dist)**2

############################################################
# Converting to SI
flux=flux*uM/(uL*uT**3)
lrange=lrange*uL

############################################################
# Converting wavelenght to micron and flux to W/m^2/micron
flux=flux*1e-6
lrange=lrange*1e6

############################################################
# Plot lambda*F_lambda vs wavelenght
plt.plot(lrange,lrange*flux,'.')
plt.xscale('log')
plt.yscale('log')
plt.show()





























































