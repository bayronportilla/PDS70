import numpy as np
import matplotlib.pyplot as plt
import sys

xval=np.linspace(0.04,120,500)
mu=0.0
sigma0=0.0059 # Keppler
Rc=40.0 # Keppler 
hp=0.1 # assian team
mb=4.0*(1.898e27/1.989e30) 
mc=4.0*(1.898e27/1.989e30) 
ms=0.76 
Nplanets=2
alpha=1e-3 # assian team
a_b=20
a_c=30.8

density=[]
randomList = ['a', 0, 2]

def udens(r):
    value=sigma0*Rc/r*np.exp(-r/Rc)
    return value

def pdens(r,mp):
    K=(mp/ms)**2*hp**-5*alpha**-1
    value=udens(r)/(1+0.04*K)
    return value

def gap_width(ap,mp):
    value=0.41*ap*(mp/ms)**0.5*hp**-0.75*alpha**-0.25
    return value

"""
# Dust Density Profile
def ddp(r,a_p):
    width=gap_width(a_p)
    wl=a_p-0.5*width
    wr=a_p+0.5*width
    if wl<=r<=wr:
        value=pdens(r)
    else:
        value=udens(r)
    return value
"""
density_profile=[]

width_b=gap_width(a_b,m)
width_c=gap_width(a_c)
wl_b=a_b-0.5*width_b
wr_b=a_b+0.5*width_b
wl_c=a_c-0.5*width_c
wr_c=a_c+0.5*width_c
for r in xval:
    if wl_b<=r<=wr_b:
        value=pdens(r)
    elif wl_c<=r<=wr_c:
        value=pdens(r)
    else:
        value=udens(r)
    density_profile.append(value)
    
    
    
density_profile=np.array(density_profile)
plt.plot(xval,density_profile)
plt.xscale('log')
plt.yscale('log')
plt.show()





    

        

sys.exit()


fig,((ax_1),(ax_2),(ax_3))=plt.subplots(1,3,figsize=(14,8))
ax_1.set_xscale('log')
ax_1.set_yscale('log')
ax_1.plot(xval,udens(xval))
ax_1.plot(xval,pdens(xval))


plt.show()
#ax_2.plot(xval,dbefore(xval))
#ax_2.plot(xval,dafter(xval,mu,sigma))
sys.exit()

#ax_2.set_xscale('log')
ax_2.set_yscale('log')
#ax_3.set_xscale('log')
ax_3.set_yscale('log')
ax_1.set_xlim(-1,1)
ax_2.set_xlim(-1,2)
ax_3.set_xlim(-1,2)
plt.show()
