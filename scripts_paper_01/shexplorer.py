import numpy as np 
import matplotlib.pyplot as plt
import sys

z1_in=0.04
z2_in=21.0
z3_in=54.0
z1_out=z2_in
z2_out=z3_in
z3_out=120.0


Npoints=1000
r_array=np.linspace(z1_in,z3_out,Npoints)
r0=100

zparams={'H0_z1':13.8,'H0_z2':13.8,'H0_z3':13.8,'b_z1':1.18,'b_z2':1.21,'b_z3':2.2}

def sh(r):
    if z1_in<=r<z1_out:
        value=zparams['H0_z1']*(r/r0)**zparams['b_z1']
    elif z2_in<=r<z2_out:
        value=zparams['H0_z2']*(r/r0)**zparams['b_z2']
    else:
        value=zparams['H0_z3']*(r/r0)**zparams['b_z3']
    return value

H_array_plus=np.array([sh(i) for i in r_array])
H_array_minus=-H_array_plus

plt.plot(r_array,H_array_plus,color='grey')
plt.plot(r_array,H_array_minus,color='grey')
plt.axvline(z1_out,linestyle='--',color='blue',linewidth=0.3)
plt.axvline(z2_out,linestyle='--',color='blue',linewidth=0.3)
plt.axvline(z3_out,linestyle='--',color='blue',linewidth=0.3)
plt.show()

    
