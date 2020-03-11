import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import fnmatch

Nzones=3

############################################################
# Working dir
folder='/data/users/bportilla/runs/final_runs/run_148/Particles'


############################################################
# Creating lists for both cases
case1=[]
case2=[]
for filename in os.listdir(folder):
    if fnmatch.fnmatch(filename,'*_f0.71_f0.29.fits.gz'):
        case1.append(filename)
    elif fnmatch.fnmatch(filename,'*_f0.58_f0.42.fits.gz'):
        case2.append(filename)


############################################################
# Define class 'Archivo'
class Archivo:
    def __init__(self,filename):
        self.filename=filename

    def zone(self):
        return self.filename[8:12]

    def binnum(self):
        return self.filename[13:17]
        
    def getEntry(self,i,j):
        hdulist=fits.open("/data/users/bportilla/runs/final_runs/run_148/Particles/%s"%self.filename)
        data=np.transpose(hdulist[0].data)
        value=data[i][j]
        return value
    
    def getWavelength(self):
        hdulist=fits.open("/data/users/bportilla/runs/final_runs/run_148/Particles/%s"%self.filename)
        data=np.transpose(hdulist[0].data)
        value=np.reshape(data[:,0:1],data.shape[0])
        return value
        

############################################################
# Create wavelength array
p1=Archivo(case1[0])
p2=Archivo(case2[0])

wl1=p1.getWavelength()
wl2=p2.getWavelength()


############################################################
# Getting Nbin for each zone
Nbin_c1_z1=0
Nbin_c1_z2=0
Nbin_c1_z3=0

Nbin_c2_z1=0
Nbin_c2_z2=0
Nbin_c2_z3=0

for particle in case1:
    p=Archivo(particle)
    if p.zone()=='0001':
        Nbin_c1_z1+=1
    elif p.zone()=='0002':
        Nbin_c1_z2+=1
    elif p.zone()=='0003':
        Nbin_c1_z3+=1

for particle in case2:
    p=Archivo(particle)
    if p.zone()=='0001':
        Nbin_c2_z1+=1
    elif p.zone()=='0002':
        Nbin_c2_z2+=1
    elif p.zone()=='0003':
        Nbin_c2_z3+=1
        
############################################################
# Create matrices to store values of Ai,Bi and Ci for case 1
# and case 2
Ac1_z1=np.zeros((len(wl1),Nbin_c1_z1))
Ac1_z2=np.zeros((len(wl1),Nbin_c1_z2))
Ac1_z3=np.zeros((len(wl1),Nbin_c1_z3))
Bc1_z1=np.zeros((len(wl1),Nbin_c1_z1))
Bc1_z2=np.zeros((len(wl1),Nbin_c1_z2))
Bc1_z3=np.zeros((len(wl1),Nbin_c1_z3))
Cc1_z1=np.zeros((len(wl1),Nbin_c1_z1))
Cc1_z2=np.zeros((len(wl1),Nbin_c1_z2))
Cc1_z3=np.zeros((len(wl1),Nbin_c1_z3))

Ac2_z1=np.zeros((len(wl2),Nbin_c2_z1))
Ac2_z2=np.zeros((len(wl2),Nbin_c2_z2))
Ac2_z3=np.zeros((len(wl2),Nbin_c2_z3))
Bc2_z1=np.zeros((len(wl2),Nbin_c2_z1))
Bc2_z2=np.zeros((len(wl2),Nbin_c2_z2))
Bc2_z3=np.zeros((len(wl2),Nbin_c2_z3))
Cc2_z1=np.zeros((len(wl2),Nbin_c2_z1))
Cc2_z2=np.zeros((len(wl2),Nbin_c2_z2))
Cc2_z3=np.zeros((len(wl2),Nbin_c2_z3))


############################################################
# Filling matrices
Nbin_index_z1=0
Nbin_index_z2=0
Nbin_index_z3=0
for particle in case1:
    p=Archivo(particle)
    if p.zone()=='0001':
        for i in range(0,len(wl1)):
            Ac1_z1[i][Nbin_index_z1]=p.getEntry(i,1)
            Bc1_z1[i][Nbin_index_z1]=p.getEntry(i,2)
            Cc1_z1[i][Nbin_index_z1]=p.getEntry(i,3)
        Nbin_index_z1+=1
    elif p.zone()=='0002':
        for i in range(0,len(wl1)):
            Ac1_z2[i][Nbin_index_z2]=p.getEntry(i,1)
            Bc1_z2[i][Nbin_index_z2]=p.getEntry(i,2)
            Cc1_z2[i][Nbin_index_z2]=p.getEntry(i,3)
        Nbin_index_z2+=1
    elif p.zone()=='0003':
        for i in range(0,len(wl1)):
            Ac1_z3[i][Nbin_index_z3]=p.getEntry(i,1)
            Bc1_z3[i][Nbin_index_z3]=p.getEntry(i,2)
            Cc1_z3[i][Nbin_index_z3]=p.getEntry(i,3)
        Nbin_index_z3+=1

#Nbin_index_z1=0
#Nbin_index_z2=0
#Nbin_index_z3=0
for particle in case2:
    p=Archivo(particle)
    if p.zone()=='0001':
        for i in range(0,len(wl2)):
            Ac2_z1[i][Nbin_index_z1]=p.getEntry(i,1)
            Bc2_z1[i][Nbin_index_z1]=p.getEntry(i,2)
            Cc2_z1[i][Nbin_index_z1]=p.getEntry(i,3)
        Nbin_index_z1+=1
    elif p.zone()=='0002':
        for i in range(0,len(wl2)):
            Ac2_z2[i][Nbin_index_z2]=p.getEntry(i,1)
            Bc2_z2[i][Nbin_index_z2]=p.getEntry(i,2)
            Cc2_z2[i][Nbin_index_z2]=p.getEntry(i,3)
        Nbin_index_z2+=1
    elif p.zone()=='0003':
        for i in range(0,len(wl2)):
            Ac2_z3[i][Nbin_index_z3]=p.getEntry(i,1)
            Bc2_z3[i][Nbin_index_z3]=p.getEntry(i,2)
            Cc2_z3[i][Nbin_index_z3]=p.getEntry(i,3)
        Nbin_index_z3+=1


############################################################
# Adding entries
sum_Ac1=np.zeros((3,len(wl1)))
sum_Bc1=np.zeros((3,len(wl1)))
sum_Cc1=np.zeros((3,len(wl1)))

sum_Ac2=np.zeros((3,len(wl1)))
sum_Bc2=np.zeros((3,len(wl1)))
sum_Cc2=np.zeros((3,len(wl1)))


for i in range(0,len(wl1)):
    sum_Ac1[0][i]=np.sum(Ac1_z1[i])
    sum_Bc1[0][i]=np.sum(Bc1_z1[i])
    sum_Cc1[0][i]=np.sum(Cc1_z1[i])

    sum_Ac1[1][i]=np.sum(Ac1_z2[i])
    sum_Bc1[1][i]=np.sum(Bc1_z2[i])
    sum_Cc1[1][i]=np.sum(Cc1_z2[i])

    sum_Ac1[2][i]=np.sum(Ac1_z3[i])
    sum_Bc1[2][i]=np.sum(Bc1_z3[i])
    sum_Cc1[2][i]=np.sum(Cc1_z3[i])
    
for i in range(0,len(wl2)):
    sum_Ac2[0][i]=np.sum(Ac2_z1[i])
    sum_Bc2[0][i]=np.sum(Bc2_z1[i])
    sum_Cc2[0][i]=np.sum(Cc2_z1[i])

    sum_Ac2[1][i]=np.sum(Ac2_z2[i])
    sum_Bc2[1][i]=np.sum(Bc2_z2[i])
    sum_Cc2[1][i]=np.sum(Cc2_z2[i])

    sum_Ac2[2][i]=np.sum(Ac2_z3[i])
    sum_Bc2[2][i]=np.sum(Bc2_z3[i])
    sum_Cc2[2][i]=np.sum(Cc2_z3[i])
    


sys.exit()
############################################################
# Plotting
fig=plt.figure(figsize=(14,12))
gs=gridspec.GridSpec(3,3)
ax1=plt.subplot(gs[0,0:1])
ax2=plt.subplot(gs[0,1:2])
ax3=plt.subplot(gs[0,2:3])
ax4=plt.subplot(gs[1,0:1])
ax5=plt.subplot(gs[1,1:2])
ax6=plt.subplot(gs[1,2:3])
ax7=plt.subplot(gs[2,0:1])
ax8=plt.subplot(gs[2,1:2])
ax9=plt.subplot(gs[2,2:3])
ax1.plot(wl1,sum_Ac1[0])
ax2.plot(wl1,sum_Bc1[0])
ax3.plot(wl1,sum_Cc1[0])
ax4.plot(wl1,sum_Ac1[1])
ax5.plot(wl1,sum_Bc1[1])
ax6.plot(wl1,sum_Cc1[1])
ax7.plot(wl1,sum_Ac1[2])
ax8.plot(wl1,sum_Bc1[2])
ax9.plot(wl1,sum_Cc1[2])
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

plt.show()



            
sys.exit()
