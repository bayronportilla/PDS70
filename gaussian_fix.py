import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.gridspec as gridspec
from scipy.interpolate import CubicSpline

def fgauss(x,sigma,mu,h):
    value=np.exp(-(x-mu)**2/(2*sigma**2))
    #return 1+h*value
    return h*(1-value)
    
xin=0.04
xcut=21.0
sigma=5
h=10
data=np.loadtxt("/data/users/bportilla/runs/run_alpha/alpha_0.001/surface_density_PDS70.dat")
r_inner=np.array([float(i) for i in data[:,0:1] if i<=xcut])



delta=np.array([2,3,5,10])
Nexperiments=len(delta)
density_matrix=np.zeros((Nexperiments,data.shape[0]))

for j in range(0,Nexperiments):
    for i in range(0,data.shape[0]):
        if i<len(r_inner):
            density_matrix[j][i]=data[i][1]*(1+fgauss(r_inner[i],sigma,xcut,delta[j]))
        else:
            density_matrix[j][i]=data[i][1]
            

fig=plt.figure(figsize=(14,8))
gs=gridspec.GridSpec(1,2)
ax1=plt.subplot(gs[0,0:1])
ax2=plt.subplot(gs[0,1:2])
for j in range(0,len(delta)):
    ax1.plot(r_inner,fgauss(r_inner,sigma,xcut,delta[j]),label='h=%.1f'%delta[j])
    ax2.plot(data[:,0:1],density_matrix[j],'--',label='h=%.1f'%delta[j])
    ax2.set_xscale('log')
    ax2.set_yscale('log')

    #ax1.set_ylim(0,1)
    #ax2.set_ylim(1e-5,1)
ax2.plot(data[:,0:1],data[:,1:2],label='initial surf. density')
ax1.legend()
ax2.legend()
plt.show()


file=open("surface_density_PDS70.dat",'w')
for i in range(0,data.shape[0]):
    file.write("%.15e %.15e\n"%(data[i][0],density_matrix[3][i]))
file.close()
