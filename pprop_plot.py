import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import fnmatch
plt.style.use('fancy')

hdulist_c1_A=fits.open('case1_A.fits') # row: wl, col: z1,z2,z3
hdulist_c1_B=fits.open('case1_B.fits')
hdulist_c1_C=fits.open('case1_C.fits')
col_A_c1=hdulist_c1_A[0].data
col_B_c1=hdulist_c1_B[0].data
col_C_c1=hdulist_c1_C[0].data

hdulist_c2_A=fits.open('case2_A.fits') # row: wl, col: z1,z2,z3
hdulist_c2_B=fits.open('case2_B.fits')
hdulist_c2_C=fits.open('case2_C.fits')
col_A_c2=hdulist_c2_A[0].data
col_B_c2=hdulist_c2_B[0].data
col_C_c2=hdulist_c2_C[0].data



############################################################
# Plotting
fig=plt.figure(figsize=(14,12))
gs=gridspec.GridSpec(3,3,hspace=0.0)
ax1=plt.subplot(gs[0,0:1])
ax2=plt.subplot(gs[0,1:2])
ax3=plt.subplot(gs[0,2:3])
ax4=plt.subplot(gs[1,0:1])
ax5=plt.subplot(gs[1,1:2])
ax6=plt.subplot(gs[1,2:3])
ax7=plt.subplot(gs[2,0:1])
ax8=plt.subplot(gs[2,1:2])
ax9=plt.subplot(gs[2,2:3])
ax1.plot(col_A_c1[:,0:1],col_A_c1[:,1:2],label='fSI=0.71, fC=0.29')
ax1.plot(col_A_c2[:,0:1],col_A_c2[:,1:2],label='fSI=0.58, fC=0.42')
ax2.plot(col_B_c1[:,0:1],col_B_c1[:,1:2])
ax2.plot(col_B_c2[:,0:1],col_B_c2[:,1:2])
ax3.plot(col_C_c1[:,0:1],col_C_c1[:,1:2])
ax3.plot(col_C_c2[:,0:1],col_C_c2[:,1:2])
ax4.plot(col_A_c1[:,0:1],col_A_c1[:,2:3])
ax4.plot(col_A_c2[:,0:1],col_A_c2[:,2:3])
ax5.plot(col_B_c1[:,0:1],col_B_c1[:,2:3])
ax5.plot(col_B_c2[:,0:1],col_B_c2[:,2:3])
ax6.plot(col_C_c1[:,0:1],col_C_c1[:,2:3])
ax6.plot(col_C_c2[:,0:1],col_C_c2[:,2:3])
ax7.plot(col_A_c1[:,0:1],col_A_c1[:,3:4])
ax7.plot(col_A_c2[:,0:1],col_A_c2[:,3:4])
ax8.plot(col_B_c1[:,0:1],col_B_c1[:,3:4])
ax8.plot(col_B_c2[:,0:1],col_B_c2[:,3:4])
ax9.plot(col_B_c1[:,0:1],col_B_c1[:,3:4])
ax9.plot(col_B_c2[:,0:1],col_B_c2[:,3:4])
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


ax7.set_xlabel(r"$\lambda \, (\mu m)$")
ax8.set_xlabel(r"$\lambda \, (\mu m)$")
ax9.set_xlabel(r"$\lambda \, (\mu m)$")

ax1.text(0.1,1e2,'zone 1, col 1')
ax4.text(0.1,1e2,'zone 2, col 1')
ax7.text(0.1,6e2,'zone 3, col 1')

ax2.text(0.1,1e2,'zone 1, col 2')
#ax5.text(0.1,1e2,'zone 2, col 2')
#ax8.text(0.1,1e2,'zone 3, col 2')

ax3.text(0.1,1e-5,'zone 1, col 3')
#ax6.text(0.1,1e0,'zone 2, col 3')
#ax9.text(0.1,1e2,'zone 3, col 3')

ax1.legend()

ax1.set_ylabel(r'$$\sum_{j=1}^{40} (\mathrm{column \, entry})_j$$')
ax4.set_ylabel(r'$$\sum_{j=1}^{40} (\mathrm{column \, entry})_j$$')
ax7.set_ylabel(r'$$\sum_{j=1}^{40} (\mathrm{column \, entry})_j$$')

fig.suptitle('opacities: $a_{\mathrm{min}}=0.001 \, (\mu m)$ and $a_{\mathrm{max}}=1000 \, (\mu m)$')
fig.savefig('opacities_amin_0.001_amax3_1000.png')
plt.show()
