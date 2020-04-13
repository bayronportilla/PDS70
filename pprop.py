import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import fnmatch

Nzones=3
Nbins=np.zeros(Nzones)
apows=np.zeros(Nzones)
psizes_min=np.zeros(Nzones)
psizes_max=np.zeros(Nzones)
psizes=[]
############################################################
# Working dir
directory="/data/users/bportilla/runs/final_runs/run_153/"
folder=directory+'Particles'
path_to_input=directory+"input.dat"
infile=open(path_to_input).readlines()
for i in range(1,Nzones+1):
    for line in infile:
        if line.split("=")[0]==("computepart0%d:ngrains"%(i)):
            Nbins[i-1]=int(line.split("=")[1])
        if line.split("=")[0]==("computepart0%d:amin"%(i)):
            psizes_min[i-1]=float(line.split("=")[1])
        if line.split("=")[0]==("computepart0%d:amax"%(i)):
            psizes_max[i-1]=float(line.split("=")[1])
        if line.split("=")[0]==("computepart0%d:apow"%(i)):
            apows[i-1]=float(line.split("=")[1])
            
for i in range(0,Nzones):
    psizes.append((psizes_min[i],psizes_max[i]))


############################################################
# Creating lists for both cases
case1=[]
case2=[]
for filename in os.listdir(folder):
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

    def f(self):
        hdulist=fits.open("/data/users/bportilla/runs/final_runs/run_153/Particles/%s"%self.filename)
        hdu=hdulist[0]
        hdr=hdu.header
        amin_bin=hdr["R_MIN"]
        amax_bin=hdr["R_MAX"]
        apow=hdr["R_POW"]
        z=int(self.zone())
        bins=Nbins[z-1]
        amin=psizes[z-1][0]
        amax=psizes[z-1][1]
        f_num=amax_bin**(-apow+4)-amin_bin**(-apow+4)
        f_den=amax**(-apow+4)-amin**(-apow+4)
        value=f_num/f_den
        return value
        
    def getEntry(self,i,j):
        hdulist=fits.open("/data/users/bportilla/runs/final_runs/run_153/Particles/%s"%self.filename)
        data=np.transpose(hdulist[0].data)
        value=data[i][j]
        return value
    
    def getWavelength(self):
        hdulist=fits.open("/data/users/bportilla/runs/final_runs/run_153/Particles/%s"%self.filename)
        data=np.transpose(hdulist[0].data)
        value=np.reshape(data[:,0:1],data.shape[0])
        return value


############################################################
# Create wavelength array
p2=Archivo(case2[0])
wl2=p2.getWavelength()


############################################################
# Getting Nbin for each zone
Nbin_c2_z1=0
Nbin_c2_z2=0
Nbin_c2_z3=0

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
for particle in case2:
    p=Archivo(particle)
    if p.zone()=='0001':
        for i in range(0,len(wl2)):
            Ac2_z1[i][Nbin_index_z1]=p.f()*p.getEntry(i,1)
            Bc2_z1[i][Nbin_index_z1]=p.f()*p.getEntry(i,2)
            Cc2_z1[i][Nbin_index_z1]=p.f()*p.getEntry(i,3)
        Nbin_index_z1+=1
    elif p.zone()=='0002':
        for i in range(0,len(wl2)):
            Ac2_z2[i][Nbin_index_z2]=p.f()*p.getEntry(i,1)
            Bc2_z2[i][Nbin_index_z2]=p.f()*p.getEntry(i,2)
            Cc2_z2[i][Nbin_index_z2]=p.f()*p.getEntry(i,3)
        Nbin_index_z2+=1
    elif p.zone()=='0003':
        for i in range(0,len(wl2)):
            Ac2_z3[i][Nbin_index_z3]=p.f()*p.getEntry(i,1)
            Bc2_z3[i][Nbin_index_z3]=p.f()*p.getEntry(i,2)
            Cc2_z3[i][Nbin_index_z3]=p.f()*p.getEntry(i,3)
        Nbin_index_z3+=1


############################################################
# Adding entries
sum_Ac2=np.zeros((3,len(wl2)))
sum_Bc2=np.zeros((3,len(wl2)))
sum_Cc2=np.zeros((3,len(wl2)))


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
# Creating HDU's for case 2
hdu_c2_A=np.stack((wl2,sum_Ac2[0],sum_Ac2[1],sum_Ac2[2]),axis=-1)
hdu_c2_B=np.stack((wl2,sum_Bc2[0],sum_Bc2[1],sum_Bc2[2]),axis=-1)
hdu_c2_C=np.stack((wl2,sum_Cc2[0],sum_Cc2[1],sum_Cc2[2]),axis=-1)
hdu_c2_A=fits.PrimaryHDU(hdu_c2_A)
hdu_c2_B=fits.PrimaryHDU(hdu_c2_B)
hdu_c2_C=fits.PrimaryHDU(hdu_c2_C)
hdu_c2_A.writeto('ext.fits')
hdu_c2_B.writeto('abs.fits')
hdu_c2_C.writeto('sca.fits')




