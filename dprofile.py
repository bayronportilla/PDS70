import numpy as np
import matplotlib.pyplot as plt
import sys

xval=np.linspace(0.04,120,500)
mu=0.0
sigma0=0.0059 # Keppler
Rc=40.0 # Keppler 
hp=0.1 # assian team
m_b=4.0*(1.898e27/1.989e30) 
m_c=4.0*(1.898e27/1.989e30) 
ms=0.76 
Nplanets=2
alpha=1e-3 # assian team
a_b=20.5
a_c=35.8

density=[]
randomList = ['a', 0, 2]

def udens(r):
    value=sigma0*Rc/r*np.exp(-r/Rc)
    return value

# perturbed density (gap_height)
def pdens(r,mp):
    K=(mp/ms)**2*hp**-5*alpha**-1
    value=udens(r)/(1+0.04*K)
    return value

"""
def gap_width(ap,mp):
    value=0.41*ap*(mp/ms)**0.5*hp**-0.75*alpha**-0.25
    return value
"""

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
r_array=[]

class Source:
    def __init__(self,mass,pos):
        self.mass=mass
        self.pos=pos
    def gap_width(self):
        value=0.41*self.pos*(self.mass/ms)**0.5*hp**-0.75*alpha**-0.25
        return value 
    def gap_height(self):
        K=(self.mass/ms)**2*hp**-5*alpha**-1
        value=udens(self.pos)/(1+0.04*K)
        return value
    def gap_limits(self):
        width=0.41*self.pos*(self.mass/ms)**0.5*hp**-0.75*alpha**-0.25
        ll=self.pos-0.5*width
        rl=self.pos+0.5*width
        return [ll,rl]

source_b=Source(m_b,a_b)
source_c=Source(m_c,a_c)

# If(source_c.gap_limits()[0]<source_b.gap_limits()[1]):
if(source_c.gap_limits()[0]<source_b.gap_limits()[1]): # Simple overlapping condition
    for r in xval:
        # Zone 1
        if r<source_b.gap_limits()[0] or r>source_c.gap_limits()[1]:
            density_profile.append(udens(r))
            r_array.append(r)
        # Zone II
        elif source_b.gap_limits()[0]<r<source_c.gap_limits()[0]:
            density_profile.append(pdens(r,source_b.mass))
            r_array.append(r)
        # Zone III
        elif source_c.gap_limits()[0]<r<source_b.gap_limits()[1]:
            value=pdens(r,source_b.mass)+pdens(r,source_c.mass)-udens(r)
            if value<0.0:
                density_profile.append(0.0)
            else:
                density_profile.append(value)
            r_array.append(r)
        # Zone IV
        elif source_b.gap_limits()[1]<r<source_c.gap_limits()[1]:
            density_profile.append(pdens(r,source_c.mass))
            r_array.append(r)
else: # No overlapping
    for r in xval:
        if source_b.gap_limits()[0]<r<source_b.gap_limits()[1]:
            density_profile.append(pdens(r,source_b.mass))
            r_array.append(r)
        elif source_c.gap_limits()[0]<r<source_c.gap_limits()[1]:
            density_profile.append(pdens(r,source_c.mass))
            r_array.append(r)
        else:
            density_profile.append(udens(r))
            r_array.append(r)
            
        
density_profile=np.array(density_profile)
r_array=np.array(r_array)

plt.plot(r_array,density_profile)
plt.xscale('log')
plt.yscale('log')
plt.ylim(-10.0,10)
plt.show()

sys.exit()

width_b=gap_width(a_b,mb)
width_c=gap_width(a_c,mc)
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
