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
folder='/data/users/bportilla/runs/final_runs/run_150/Particles'


############################################################
# Creating lists for both cases
case1=[]
case2=[]
for filename in os.listdir(folder):
    #if fnmatch.fnmatch(filename,'*_f0.71_f0.29.fits.gz'):
    #    case1.append(filename)
    if fnmatch.fnmatch(filename,'*_f0.58_f0.42.fits.gz'):
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
        hdulist=fits.open("/data/users/bportilla/runs/final_runs/run_150/Particles/%s"%self.filename)
        data=np.transpose(hdulist[0].data)
        value=data[i][j]
        return value
    
    def getWavelength(self):
        hdulist=fits.open("/data/users/bportilla/runs/final_runs/run_150/Particles/%s"%self.filename)
        data=np.transpose(hdulist[0].data)
        value=np.reshape(data[:,0:1],data.shape[0])
        return value
        

############################################################
# Create wavelength array
#p1=Archivo(case1[0])
p2=Archivo(case2[0])

#wl1=p1.getWavelength()
wl2=p2.getWavelength()


############################################################
# Getting Nbin for each zone
"""
Nbin_c1_z1=0
Nbin_c1_z2=0
Nbin_c1_z3=0
"""
Nbin_c2_z1=0
Nbin_c2_z2=0
Nbin_c2_z3=0

"""
for particle in case1:
    p=Archivo(particle)
    if p.zone()=='0001':
        Nbin_c1_z1+=1
    elif p.zone()=='0002':
        Nbin_c1_z2+=1
    elif p.zone()=='0003':
        Nbin_c1_z3+=1
"""
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
"""
Ac1_z1=np.zeros((len(wl1),Nbin_c1_z1))
Ac1_z2=np.zeros((len(wl1),Nbin_c1_z2))
Ac1_z3=np.zeros((len(wl1),Nbin_c1_z3))
Bc1_z1=np.zeros((len(wl1),Nbin_c1_z1))
Bc1_z2=np.zeros((len(wl1),Nbin_c1_z2))
Bc1_z3=np.zeros((len(wl1),Nbin_c1_z3))
Cc1_z1=np.zeros((len(wl1),Nbin_c1_z1))
Cc1_z2=np.zeros((len(wl1),Nbin_c1_z2))
Cc1_z3=np.zeros((len(wl1),Nbin_c1_z3))
"""
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
"""
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
"""
Nbin_index_z1=0
Nbin_index_z2=0
Nbin_index_z3=0
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
"""
sum_Ac1=np.zeros((3,len(wl1)))
sum_Bc1=np.zeros((3,len(wl1)))
sum_Cc1=np.zeros((3,len(wl1)))
"""
sum_Ac2=np.zeros((3,len(wl2)))
sum_Bc2=np.zeros((3,len(wl2)))
sum_Cc2=np.zeros((3,len(wl2)))

"""
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
"""

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


############################################################
# Creating HDU's for case 1 and case 2
"""
hdu_c1_A=np.stack((wl1,sum_Ac1[0],sum_Ac1[1],sum_Ac1[2]),axis=-1)
hdu_c1_B=np.stack((wl1,sum_Bc1[0],sum_Bc1[1],sum_Bc1[2]),axis=-1)
hdu_c1_C=np.stack((wl1,sum_Cc1[0],sum_Cc1[1],sum_Cc1[2]),axis=-1)
hdu_c1_A=fits.PrimaryHDU(hdu_c1_A)
hdu_c1_B=fits.PrimaryHDU(hdu_c1_B)
hdu_c1_C=fits.PrimaryHDU(hdu_c1_C)
hdu_c1_A.writeto('case1_A.fits')
hdu_c1_B.writeto('case1_B.fits')
hdu_c1_C.writeto('case1_C.fits')
"""
hdu_c2_A=np.stack((wl2,sum_Ac2[0],sum_Ac2[1],sum_Ac2[2]),axis=-1)
hdu_c2_B=np.stack((wl2,sum_Bc2[0],sum_Bc2[1],sum_Bc2[2]),axis=-1)
hdu_c2_C=np.stack((wl2,sum_Cc2[0],sum_Cc2[1],sum_Cc2[2]),axis=-1)
hdu_c2_A=fits.PrimaryHDU(hdu_c2_A)
hdu_c2_B=fits.PrimaryHDU(hdu_c2_B)
hdu_c2_C=fits.PrimaryHDU(hdu_c2_C)
hdu_c2_A.writeto('case2_A.fits')
hdu_c2_B.writeto('case2_B.fits')
hdu_c2_C.writeto('case2_C.fits')




