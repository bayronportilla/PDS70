from mcmax3d_analysis.mcmax3d_observables import convert_flux
import numpy as np
import matplotlib.pyplot as plt
import sys
plt.style.use('fancy')

N_sim=88
############################################################
# Converting MCMax3D spectra 
#convert_flux('/data/users/bportilla/runs/my_MCMax3D_run001/output/MCSpec0001.dat','spectrum_PDS70_system.dat')
#convert_flux('/data/users/bportilla/runs/my_MCMax3D_run001/output/star0001.dat','spectrum_PDS70_star.dat')
convert_flux('/data/users/bportilla/runs/PDS70_runs/run_%d/output/MCSpec0001.dat'%N_sim,'spectrum_PDS70_system.dat')
convert_flux('/data/users/bportilla/runs/PDS70_runs/run_%d/output/star0001.dat'%N_sim,'spectrum_PDS70_star.dat')

############################################################
# Loading surface density profile
#surface_density=np.loadtxt('/data/users/bportilla/runs/my_MCMax3D_run001/surface_density_PDS70_8.0_0.0_30.0.dat')
surface_density=np.loadtxt('/data/users/bportilla/runs/PDS70_runs/run_%d/surface_density_PDS70.dat'%N_sim)

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

############################################################
# Plot lambda*F_lambda vs wavelenght
lw=1.5
pw=7
fig,(ax_1,ax_2) = plt.subplots(1, 2,figsize=(14,5),sharex=False)
ax_1.plot(x_system,x_system*y_system,label='star+disk+cpd',linewidth=lw)
ax_1.plot(x_star,x_star*y_star,label='only star')
ax_1.errorbar(x_photometric,x_photometric*y_photometric,yerr=y_system_error,label='photometric data',fmt='.',markersize=pw)
ax_1.legend()
ax_1.set_xscale('log')
ax_1.set_yscale('log')
ax_1.set_xlabel(r'$\lambda \, (\mu\mathrm{m})$')
ax_1.set_ylabel(r'$\lambda F_{\lambda}$ (W/m^2)')
ax_1.set_ylim(1e-17,1e-12)
ax_1.set_xlim(0.1,2e3)

ax_2.plot(surface_density[:,0:1],surface_density[:,1:2])
ax_2.set_xscale('log')
ax_2.set_yscale('log')
ax_2.set_xlabel(r'$r$ (AU)')
ax_2.set_ylabel(r'$\Sigma_{\mathrm{dust}}$ (g/cm^2)')
ax_2.set_ylim(1e-6,1e3)

#fig.savefig('new/PDS70_MCMax3D_spectrum_3M_56.png')
plt.show()



