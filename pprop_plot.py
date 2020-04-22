import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import fnmatch
plt.style.use('fancy')

hdulist_ext=fits.open('ext.fits') # row: wl, col: z1,z2,z3
hdulist_abso=fits.open('abs.fits')
hdulist_sca=fits.open('sca.fits')
ext=hdulist_ext[0].data
abso=hdulist_abso[0].data
sca=hdulist_sca[0].data
Nzones=1#ext.shape[1]-1


############################################################
# Plotting
fig=plt.figure(figsize=(9,6))
gs=gridspec.GridSpec(1,Nzones,hspace=0.0)
for i in range(0,Nzones):
    ax=plt.subplot(gs[0,i])
    ax.plot(ext[:,0:1],ext[:,1:2],label=r"$\kappa_{\nu}^{\mathrm{ext}}$")
    ax.plot(abso[:,0:1],abso[:,1:2],label=r"$\kappa_{\nu}^{\mathrm{abs}}$")
    ax.plot(sca[:,0:1],sca[:,1:2],label=r"$\kappa_{\nu}^{\mathrm{sca}}$")
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
ax.legend(fontsize=18)
ax.set_xlabel(r"$\lambda \, (\mu m)$")
ax.set_ylabel(r"Dust oppacities (cm$^2$/g(dust))")
plt.show()
fig.savefig("opa3.png")
    
"""
ax1=plt.subplot(gs[0,0:1])
ax2=plt.subplot(gs[0,1:2])
ax3=plt.subplot(gs[0,2:3])
ax4=plt.subplot(gs[1,0:1])
ax5=plt.subplot(gs[1,1:2])
ax6=plt.subplot(gs[1,2:3])
ax7=plt.subplot(gs[2,0:1])
ax8=plt.subplot(gs[2,1:2])
ax9=plt.subplot(gs[2,2:3])
ax10=plt.subplot(gs[3,0:1])
ax11=plt.subplot(gs[3,1:2])
ax12=plt.subplot(gs[3,2:3])
"""
"""
ax1.plot(ext[:,0:1],ext[:,1:2])#,label='fSI=0.58, fC=0.42')
ax2.plot(ext[:,0:1],ext[:,2:3])
ax3.plot(ext[:,0:1],ext[:,3:4])
ax4.plot(abso[:,0:1],abso[:,1:2])
ax5.plot(abso[:,0:1],abso[:,2:3])
ax6.plot(abso[:,0:1],abso[:,3:4])
ax7.plot(sca[:,0:1],sca[:,1:2])
ax8.plot(sca[:,0:1],sca[:,2:3])
ax9.plot(sca[:,0:1],sca[:,3:4])
ax10.plot(ext[:,0:1],ext[:,1:2]-(abso[:,1:2]+sca[:,1:2]))
ax11.plot(ext[:,0:1],ext[:,2:3]-(abso[:,2:3]+sca[:,2:3]))
ax12.plot(ext[:,0:1],ext[:,3:4]-(abso[:,3:4]+sca[:,3:4]))
"""

sys.exit()

############################################################
# Plotting
fig=plt.figure(figsize=(14,12))
gs=gridspec.GridSpec(4,3,hspace=0.0)
ax1=plt.subplot(gs[0,0])
ax2=plt.subplot(gs[0,1])
ax3=plt.subplot(gs[0,2])
ax4=plt.subplot(gs[1,0])
ax5=plt.subplot(gs[1,1])
ax6=plt.subplot(gs[1,2])
ax7=plt.subplot(gs[2,0])
ax8=plt.subplot(gs[2,1])
ax9=plt.subplot(gs[2,2])
ax10=plt.subplot(gs[3,0])
ax11=plt.subplot(gs[3,1])
ax12=plt.subplot(gs[3,2])

ax1.plot(ext[:,0:1],ext[:,1:2])#,label='fSI=0.58, fC=0.42')
ax2.plot(ext[:,0:1],ext[:,2:3])
ax3.plot(ext[:,0:1],ext[:,3:4])
ax4.plot(abso[:,0:1],abso[:,1:2])
ax5.plot(abso[:,0:1],abso[:,2:3])
ax6.plot(abso[:,0:1],abso[:,3:4])
ax7.plot(sca[:,0:1],sca[:,1:2])
ax8.plot(sca[:,0:1],sca[:,2:3])
ax9.plot(sca[:,0:1],sca[:,3:4])
ax10.plot(ext[:,0:1],ext[:,1:2]-(abso[:,1:2]+sca[:,1:2]))
ax11.plot(ext[:,0:1],ext[:,2:3]-(abso[:,2:3]+sca[:,2:3]))
ax12.plot(ext[:,0:1],ext[:,3:4]-(abso[:,3:4]+sca[:,3:4]))


ax1.set_xscale('log')
ax2.set_xscale('log')
ax3.set_xscale('log')
ax4.set_xscale('log')
ax5.set_xscale('log')
ax6.set_xscale('log')
ax7.set_xscale('log')
ax8.set_xscale('log')
ax9.set_xscale('log')
ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
ax4.set_yscale('log')
ax5.set_yscale('log')
ax6.set_yscale('log')
ax7.set_yscale('log')
ax8.set_yscale('log')
ax9.set_yscale('log')

ax10.set_xlabel(r"$\lambda \, (\mu m)$")
ax11.set_xlabel(r"$\lambda \, (\mu m)$")
ax12.set_xlabel(r"$\lambda \, (\mu m)$")

ax1.set_title('zone 1')
ax2.set_title('zone 2')
ax3.set_title('zone 3')


ax1.legend()

ax1.set_ylabel(r'$\kappa \, (\mathrm{cm}^2/\mathrm{g})$')
ax4.set_ylabel(r'$\alpha \, (\mathrm{cm}^2/\mathrm{g})$')
ax7.set_ylabel(r'$\sigma \, (\mathrm{cm}^2/\mathrm{g})$')

#fig.suptitle('opacities: $a_{\mathrm{min},3}=0.001 \, (\mu m)$ and $a_{\mathrm{max},3}=1000 \, (\mu m)$')
fig.savefig('opacities_run_153.png')
plt.show()
