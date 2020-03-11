import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import sys
import fnmatch

Nzones=3
Nbin=40

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
Ac1_z1=np.zeros((len(wl1),Nbin_c1_z1))
Ac1_z2=np.zeros((len(wl1),Nbin_c1_z2))
Ac1_z3=np.zeros((len(wl1),Nbin_c1_z3))

Bc1_z1=np.zeros((len(wl1),Nbin_c1_z1))
Bc1_z2=np.zeros((len(wl1),Nbin_c1_z2))
Bc1_z3=np.zeros((len(wl1),Nbin_c1_z3))

Cc1_z1=np.zeros((len(wl1),Nbin_c1_z1))
Cc1_z2=np.zeros((len(wl1),Nbin_c1_z2))
Cc1_z3=np.zeros((len(wl1),Nbin_c1_z3))


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

print(Ac1_z1)
            
sys.exit()
for particle in case2:
    p=Archivo(particle)
    
    print(p.zone())

sys.exit()
f1=Archivo(case1[0])
print(f1.zone())
print(f1.getEntry(-1,0))


sys.exit()
hdulist=fits.open("/data/users/bportilla/runs/final_runs/run_148/Particles/particle0002_0006_0001_f0.71_f0.29.fits.gz")
hdulist.info()
data=hdulist[0].data
print(np.transpose(data)[:,0:1])
