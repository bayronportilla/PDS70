from mcmax3d_analysis.mcmax3d_observables import convert_flux
import numpy as np
import matplotlib.pyplot as plt
import sys
plt.style.use('fancy')

N_files=len(sys.argv)-1

for i in range(1,len(sys.argv)):

    ############################################################
    # Converting MCMax3D spectra 
    convert_flux('/data/users/bportilla/runs/PDS70_runs/run_%d/output/MCSpec0001.dat'%(int(sys.argv[i])),'spectrum_PDS70_system_%d.dat'%(int(sys.argv[i])))
    convert_flux('/data/users/bportilla/runs/PDS70_runs/run_%d/output/star0001.dat'%(int(sys.argv[i])),'spectrum_PDS70_star_%d.dat'%(int(sys.argv[i])))

    ############################################################
    # Loading data (converted) from MCMax3D 
    system_spectrum=np.loadtxt('spectrum_PDS70_system_%d.dat'%(int(sys.argv[i])))
    stellar_spectrum=np.loadtxt('spectrum_PDS70_star_%d.dat'%(int(sys.argv[i])))

    ############################################################
    # Loading surface density profile
    surface_density=np.loadtxt('/data/users/bportilla/runs/PDS70_runs/run_%d/surface_density_PDS70.dat'%(int(sys.argv[i])))

    if i==1:
        ############################################################
        # Creating arrays for plotting
        x_system=np.reshape(system_spectrum[:,0:1],system_spectrum.shape[0])
        y_system=np.reshape(system_spectrum[:,1:2],system_spectrum.shape[0])
        x_star=np.reshape(stellar_spectrum[:,0:1],stellar_spectrum.shape[0])
        y_star=np.reshape(stellar_spectrum[:,1:2],stellar_spectrum.shape[0])
    else:
        x_system=np.vstack((x_system,np.reshape(system_spectrum[:,0:1],system_spectrum.shape[0])))
        y_system=np.vstack((y_system,np.reshape(system_spectrum[:,1:2],system_spectrum.shape[0])))
        x_star=np.vstack((x_star,np.reshape(stellar_spectrum[:,0:1],stellar_spectrum.shape[0])))
        y_star=np.vstack((y_star,np.reshape(stellar_spectrum[:,1:2],stellar_spectrum.shape[0])))
    

############################################################
# Loading photometric data
photometric_data=np.loadtxt('PDS70_photometry.dat')
x_photometric=np.reshape(photometric_data[:,0:1],photometric_data.shape[0])
y_photometric=np.reshape(photometric_data[:,1:2],photometric_data.shape[0])
y_photometric_error=np.reshape(photometric_data[:,2:3],photometric_data.shape[0])

############################################################
# Creating array of errors
y_system_error=x_photometric*y_photometric_error


############################################################
# Plot lambda*F_lambda vs wavelenght
fig,(ax_1,ax_2) = plt.subplots(1, 2,figsize=(14,5),sharex=False)
for i in range(0,N_files):
    ax_1.plot(x_system[i],x_system[i]*y_system[i],label=('run %d'%(i+1)))
    ax_1.plot(x_star[i],x_star[i]*y_star[i])
ax_1.errorbar(x_photometric,x_photometric*y_photometric,yerr=y_system_error,label='photometric data',fmt='.',markersize=5.0)
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



