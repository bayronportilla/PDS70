from mcmax3d_analysis.mcmax3d_observables import convert_flux
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('fancy')

############################################################
# Converting MCMax3D spectra 
convert_flux('/data/users/bportilla/runs/my_MCMax3D_run001/output/MCSpec0001.dat','spectrum_PDS70_system.dat')
convert_flux('/data/users/bportilla/runs/my_MCMax3D_run001/output/star0001.dat','spectrum_PDS70_star.dat')

############################################################
# Loading data 
system_spectrum=np.loadtxt('spectrum_PDS70_system.dat')
stellar_spectrum=np.loadtxt('spectrum_PDS70_star.dat')

############################################################
# Creating arrays for plotting
x_system=np.reshape(system_spectrum[:,0:1],system_spectrum.shape[0])
y_system=np.reshape(system_spectrum[:,1:2],system_spectrum.shape[0])
x_star=np.reshape(stellar_spectrum[:,0:1],stellar_spectrum.shape[0])
y_star=np.reshape(stellar_spectrum[:,1:2],stellar_spectrum.shape[0])

############################################################
# Plot lambda*F_lambda vs wavelenght
fig=plt.figure(figsize=(16,8))
ax_1=plt.axes()
ax_1.plot(x_system,x_system*y_system,label='star+disk+cpd')
ax_1.plot(x_star,x_star*y_star,label='only star')
ax_1.legend()
ax_1.set_xscale('log')
ax_1.set_yscale('log')
ax_1.set_xlabel(r'$\lambda \, (\mu\mathrm{m})$')
ax_1.set_ylabel(r'$\lambda F_{\lambda}$ (W/m^2)')
ax_1.set_ylim(1e-17,1e-12)
ax_1.set_xlim(0.1,2e3)

fig.savefig('whole_system.png')
plt.show()



