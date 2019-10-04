import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy import constants as cte
from scipy.integrate import romb
from chiang_goldreich_model import *
from height_scale import *
plt.style.use('fancy')

############################################################
# Constants and definitions
m_a=0.76 # Stellar mass (solar masses)
m_b=(7.5)*cte.M_jup.value/cte.M_sun.value # Mass of PDS70 b (solar masses)
m_c=(8.0)*cte.M_jup.value/cte.M_sun.value # Mass of PDS70 c (solar masses)
a_b=20.6 # Semi-major axis of PDS70 b (AU)
a_c=34.5 # Semi-major axis of PDS70 c (AU)
mu=0.0
alpha=1e-3

class Companion:
    def __init__(self,mass,pos):
        self.mass=mass
        self.pos=pos
    def gap_width(self):
        value=0.41*self.pos*(self.mass/m_a)**0.5*(h_scale(m_a,self.pos)/self.pos)**(-3./4)*alpha**(-1/4.)
        return value 

############################################################
# Case 1. No gaps
R_in=0.04 # Disk's inner limit (AU)
R_out=120.0 # Disk's outer limit (AU)
R_exp=40.0 # Characteristic radius (AU)
N_points=1.025e3 # Number of discrete points where density is computed
M_dust=3.0e-5 # Total dust mass (solar masses)
r_array=np.linspace(R_in,R_out,N_points)

def surface_density(r):
    
    ############################################################
    #
    # Returns the unperturbed density (that with no gaps) 
    # in g/cm2. The input is the radial distance in AU.
    #
    ############################################################

    exp_value=np.exp(-R_in/R_exp)-np.exp(-R_out/R_exp)
    sigma_0=M_dust/(2.0*np.pi*R_exp**2*exp_value) * (cte.M_sun.value*1000/(cte.au.value*100)**2)
    value=sigma_0*R_exp/r*np.exp(-r/R_exp)
    return value

"""
rho_array=surface_density(r_array)

file=open('surface_density_PDS70.dat','w')
for i in range(0,len(r_array)):
    file.write('%.15e %.15e\n'%(r_array[i],rho_array[i])) # col1:  r (AU), col2: surf. density (g/cm2)

x_integrate=r_array*cte.au.value*100
y_integrate=x_integrate*rho_array
dx=x_integrate[1]-x_integrate[0]

print(2*np.pi*romb(y_integrate,dx)/(cte.M_sun.value*1000))

plt.plot(r_array,rho_array)
plt.xscale('log')
plt.yscale('log')
plt.show()


sys.exit()

############################################################
# Case 2. One gap
def surface_density_one_gap(r,m,a):
    planet_b=Companion(m,a)
    K=(m/m_a)**2*(h_scale(m_a,a)/a)**-5*alpha**-1
    d_min_b=a-0.5*planet_b.gap_width()
    d_max_b=a+0.5*planet_b.gap_width()
    exp_value= np.exp(-R_in/R_exp)-np.exp(-R_out/R_exp)+\
               np.exp(-d_max_b/R_exp)-np.exp(-d_min_b/R_exp)+\
               (np.exp(-d_min_b/R_exp)-np.exp(d_max_b/R_exp))/(1+0.04*K)
    sigma_0=M_dust/(2.0*np.pi*R_exp**2*exp_value) * (cte.M_sun.value*1000/(cte.au.value*100)**2)
    try:
        if R_in<=r<d_min_b or d_max_b<r<=R_out:
            value=sigma_0*R_exp/r*np.exp(-r/R_exp)
            return value
        elif d_min_b<r<d_max_b :
            value=(sigma_0*R_exp/r*np.exp(-r/R_exp))/(1+0.04*K)
            return value
        else:
            raise ValueError
    except ValueError:
        return('r is not between R_in and R_out')                                    

rho_array=[]
for r in r_array:
    rho_array.append(surface_density_one_gap(r,m_b,a_b))
rho_array=np.array(rho_array)

x_integrate=r_array*cte.au.value*100
y_integrate=x_integrate*rho_array
dx=x_integrate[1]-x_integrate[0]
print(2*np.pi*romb(y_integrate,dx)/(cte.M_sun.value*1000))


plt.plot(r_array,rho_array)
plt.xscale('log')
plt.yscale('log')
plt.show()
"""



############################################################
# Case 3. Two gaps overlapped 
def surface_density_two_gaps_overlapped(r,m_b,a_b,m_c,a_c):
    planet_b=Companion(m_b,a_b)
    planet_c=Companion(m_c,a_c)
    K_b=(m_b/m_a)**2*(h_scale(m_a,a_b)/a_b)**-5*alpha**-1
    K_c=(m_c/m_a)**2*(h_scale(m_a,a_c)/a_c)**-5*alpha**-1
    d_min_b=a_b-0.5*planet_b.gap_width()
    d_max_b=a_b+0.5*planet_b.gap_width()
    d_min_c=a_c-0.5*planet_c.gap_width()
    d_max_c=a_c+0.5*planet_c.gap_width()
    delta=1e-20 # Depletion factor (dimensionless)
    D_1=np.exp(-R_in/R_exp)-np.exp(-d_min_b/R_exp)
    D_2=np.exp(-d_min_b/R_exp)-np.exp(-d_min_c/R_exp)
    D_3=np.exp(-d_min_c/R_exp)-np.exp(-d_max_b/R_exp)
    D_4=np.exp(-d_max_b/R_exp)-np.exp(-d_max_c/R_exp)
    D_5=np.exp(-d_max_c/R_exp)-np.exp(-R_out/R_exp)
    exp_value=D_1+D_2/(1+0.04*K_b)+delta*D_3+D_4/(1+0.04*K_c)+D_5
    sigma_0=M_dust/(2.0*np.pi*R_exp**2*exp_value) * (cte.M_sun.value*1000/(cte.au.value*100)**2)
    #value=sigma_0*R_exp/r*np.exp(-r/R_exp)
    try:
        if R_in<=r<d_min_b or d_max_c<r<=R_out:
            value=sigma_0*R_exp/r*np.exp(-r/R_exp)
            return value
        elif d_min_b<r<d_min_c:
            value=(sigma_0*R_exp/r*np.exp(-r/R_exp))/(1+0.04*K_b)
            return value
        elif d_min_c<r<d_max_b:
            value=delta*sigma_0*R_exp/r*np.exp(-r/R_exp)
            return value
        elif d_max_b<r<d_max_c:
            value=(sigma_0*R_exp/r*np.exp(-r/R_exp))/(1+0.04*K_c)
            return value
        else:
            raise ValueError
    except ValueError:
        return('r is not between R_in and R_out')                                    
    
rho_array=[]
for r in r_array:
    rho_array.append(surface_density_two_gaps_overlapped(r,m_b,a_b,m_c,a_c))
rho_array=np.array(rho_array)

"""
x_integrate=r_array*cte.au.value*100
y_integrate=x_integrate*rho_array
dx=x_integrate[1]-x_integrate[0]
print(2*np.pi*romb(y_integrate,dx)/(cte.M_sun.value*1000))
file=open('surface_density_PDS70.dat','w')
for i in range(0,len(r_array)):
    file.write('%.15e %.15e\n'%(r_array[i],rho_array[i])) # col1:  r (AU), col2: surf. density (g/cm2)
""" 

plt.plot(r_array,rho_array)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$r$(AU)')
plt.ylabel('$\Sigma_{\mathrm{dust}}$(g/cm$^2$)')
plt.ylim(1e-6,1e3)
plt.savefig('surface_density_PDS70.png')
plt.show()


sys.exit()








sys.exit()


xval=np.linspace(0.04,120,500)
mu=0.0
sigma0=0.0059 # Keppler
R_exp=40.0 # Keppler 
hp=0.1 # assian team
m_b=7.5*(1.898e27/1.989e30) 
m_c=8.0*(1.898e27/1.989e30) 
ms=0.76 
Nplanets=2
alpha=1e-3 # assian team
a_b=20.6
a_c=34.5
beta=1.25
h100=10


sys.exit()

def h(r):
    value=h100*(r/100)**beta
    return value

density=[]
randomList = ['a', 0, 2]

def udens(r):
    value=sigma0*R_exp/r*np.exp(-r/R_exp)
    return value

# perturbed density (gap_height)
def pdens(r,mp):
    K=(mp/ms)**2*hp**-5*alpha**-1
    #K=(mp/ms)**2*h(r)**-5*alpha**-1
    value=udens(r)/(1+0.04*K)
    return value

density_profile=[]
r_array=[]

class Companion:
    def __init__(self,mass,pos):
        self.mass=mass
        self.pos=pos
    def gap_width(self):
        #value=0.41*self.pos*(self.mass/ms)**0.5*h(self.pos)**-0.75*alpha**-0.25
        value=0.41*self.pos*(self.mass/ms)**0.5*hp**-0.75*alpha**-0.25
        return value 
    def gap_height(self):
        #K=(self.mass/ms)**2*h(self.pos)**-5*alpha**-1
        K=(self.mass/ms)**2*hp**-5*alpha**-1
        value=udens(self.pos)/(1+0.04*K)
        return value
    def gap_limits(self):
        #width=0.41*self.pos*(self.mass/ms)**0.5*h(self.pos)**-0.75*alpha**-0.25
        width=0.41*self.pos*(self.mass/ms)**0.5*hp**-0.75*alpha**-0.25
        ll=self.pos-0.5*width
        rl=self.pos+0.5*width
        return [ll,rl]

source_b=Companion(m_b,a_b)
source_c=Companion(m_c,a_c)

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


############################################################
# Printing information
print('limits source b:',(source_b.gap_limits()))
print('limits source c:',(source_c.gap_limits()))

sys.exit()




density_profile=np.array(density_profile)
r_array=np.array(r_array)
            
############################################################
# Fancy style
plt.rc('font', **{'family': 'family', 'family': ['serif']})
plt.rc('text', usetex=True) #uses Latex instead of Tex to compile axes labels

print(source_b.gap_limits()[0])
#sys.exit()
fig=plt.figure()
ax=plt.axes()
ax.plot(r_array,density_profile)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$r$ (AU)')
ax.set_ylabel('$\Sigma (r)$ (g/cm$^2$)')
ax.axvspan(source_b.gap_limits()[0],source_b.gap_limits()[1], ymin=0.0, ymax=10.0, alpha=0.1, color='blue')
ax.axvspan(source_c.gap_limits()[0],source_c.gap_limits()[1], ymin=0.0, ymax=10.0, alpha=0.1, color='green')
#ax.fill_between(r_array,0,10, where=r_array<33,facecolor='green', alpha=0.5)
#ax.fill_between(r_array,0,10, where=r_array>0.8,facecolor='red', alpha=0.1)
ax.axvline(a_b,linestyle='--',linewidth=0.5,color='black')
ax.axvline(a_c,linestyle='--',linewidth=0.5,color='black')
ax.set_ylim(-10.0,10)
#ax.tick_params(direction='in')
fig.savefig('dprofile_PDS70.png')
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
