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

"""
############################################################
# Synthetic spectrum
synthetic_spec=np.loadtxt('Kurucz_04000_4.0.dat',skiprows=3)
#synthetic_spec=np.loadtxt('Kurucz_04000_0.0.dat',skiprows=3)

x_synthetic=[]
y_synthetic=[]
for i in range(0,synthetic_spec.shape[0]):
    x_synthetic.append(synthetic_spec[i][0])
    y_synthetic.append(synthetic_spec[i][1])
x_synthetic=np.array(x_synthetic)
y_synthetic=np.array(y_synthetic) # erg/cm^2/s/Hz/sr
y_synthetic=y_synthetic*np.pi # erg/cm^2/s/Hz
y_synthetic=y_synthetic/300.0e6 # erg/s*cm^2*m
y_synthetic=y_synthetic/1e3 # W/m^2*m
y_synthetic=y_synthetic*1e-6 # W/m^2*microns
x_synthetic=x_synthetic/1000.0 # microns

y_plot=x_synthetic*y_synthetic*(R/dist)**2


plt.plot(x_synthetic,y_synthetic*x_synthetic*(R/dist)**2)
plt.xscale('log')
plt.yscale('log')
plt.show()

#sys.exit()
"""
############################################################
# Wavelenght range (microns)
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
# Generating file of theoretical spectrum for MCMax3D
file=open('black_body_PDS70.dat','w')
for i in range(0,len(lrange)):
    file.write("%.16e %.16e\n"%(lrange[i],(flux*lrange)[i]))

sys.exit()
############################################################
# Plot lambda*F_lambda vs wavelenght
plt.plot(lrange,flux*lrange)
#plt.plot(lrange,flux,'.')
#plt.plot(x_synthetic,y_synthetic*x_synthetic*(R/dist)**2*1e29)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\lambda$ (microns)')
plt.ylabel('$\lambda F_{\lambda}$ (W/m^2)')
plt.xlim(1e-1,1e3)
plt.ylim(1e-17,1e-12)
plt.show()





























































