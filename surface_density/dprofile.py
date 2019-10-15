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
R_in=0.04 # Disk's inner limit (AU)
R_out=120.0 # Disk's outer limit (AU)
R_exp=40.0 # Characteristic radius (AU)
k=10
N_points=2**k+1 # Number of discrete points where density is computed. (=2^k+1, k integer)
M_dust=3.0e-4 # Total dust mass (solar masses)
g=1.0 # Ratio M_disk/M_dust
M_dust_disk=g*M_dust # Dust mass in the disk
r_array=np.linspace(R_in,R_out,N_points)


############################################################
# Define class Companion for the planets
class Companion:
    def __init__(self,mass,pos):
        self.mass=mass
        self.pos=pos
    def gap_width(self):
        value=0.41*self.pos*(self.mass/m_a)**0.5*(h_scale(m_a,self.pos)/self.pos)**(-3./4)*alpha**(-1/4.)
        return value 

############################################################
# Verification of the total dust mass retrieved by the 
# density profile.
def verification_total_mass(m,x,y):

    ############################################################
    #
    # This function performs the surface integral of the density
    # profile. Inputs are: the total dust mass (the constraint) in
    # solar masses, the array of distances (in AU) and the array
    # with the surface density profile array (in g/cm^2). The 
    # output is the error between the total dust mass generated
    # by the density profile and the actual value. Therefore, 
    # if the density profile reproduces exactly the value of 
    # the constraint, the value returned will be zero. 
    #
    # IMPORTANT: the number of points of the density profile array
    # can not be arbitrary. It should respect the rule: N=2^k+1
    # with k integer. Also, the values of the independent variable
    # i.e. the radial distance, must be equally spaced. 
    #
    ############################################################
    
    x_integrate=x*cte.au.value*100
    y_integrate=x_integrate*y
    dx=x_integrate[1]-x_integrate[0]
    dust_mass_integrated=2*np.pi*romb(y_integrate,dx)/(cte.M_sun.value*1000)
    #dust_mass_integrated=romb(y_integrate,dx)/(cte.M_sun.value*1000)
    error=abs(m-dust_mass_integrated)
    #return dust_mass_integrated
    return error


############################################################
# Case 0. Fictitious profile (one gap)
def fictitious_density_profile(r,d_min,d_max):
    delta_1=0.05
    delta_2=1e-15 # Depletion factor (dimensionless)
    delta_3=1.0
    exp_value= delta_1*(np.exp(-R_in/R_exp)-np.exp(-d_min/R_exp))+\
               delta_2*(np.exp(-d_min/R_exp)-np.exp(-d_max/R_exp))+\
               delta_3*(np.exp(-d_max/R_exp)-np.exp(-R_out/R_exp))
    
    sigma_0=M_dust_disk/(2.0*np.pi*R_exp**2*exp_value) * (cte.M_sun.value*1000/(cte.au.value*100)**2)
    #sigma_0=M_dust_disk/(R_exp**2*exp_value) * (cte.M_sun.value*1000/(cte.au.value*100)**2)
    div=1
    try:
        if R_in<=r<d_min:
            value=delta_1*sigma_0*R_exp/r*np.exp(-r/R_exp)/div
            return value
        elif d_min<r<d_max :
            value=delta_2*(sigma_0*R_exp/r*np.exp(-r/R_exp))/div
            return value
        elif d_max<=r<=R_out:
            value=delta_3*sigma_0*R_exp/r*np.exp(-r/R_exp)/div
            return value
        else:
            raise ValueError
    except ValueError:
        return('r is not between R_in and R_out')                                    

dmin=20
dmax=30
############################################################
# Generating profiles
rho_array=[]
for r in r_array:
    #rho_array.append(surface_density_two_gaps_overlapped(r,m_b,a_b,m_c,a_c))
    rho_array.append(fictitious_density_profile(r,dmin,dmax))
rho_array=np.array(rho_array)

plt.plot(r_array,rho_array)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$r$(AU)')
plt.ylabel('$\Sigma_{\mathrm{dust}}$(g/cm$^2$)')
plt.ylim(1e-6,1e3)
#plt.savefig('surface_density_PDS70.png')
plt.show()

print("integrated mass:",verification_total_mass(M_dust_disk,r_array,rho_array))
print("M dust disk:",M_dust_disk)
print(verification_total_mass(M_dust_disk,r_array,rho_array)-M_dust_disk)


############################################################
# Saving data
file=open('surface_density_PDS70_case_2.dat','w')
if verification_total_mass(M_dust_disk,r_array,rho_array)<1e-6:
    for i in range(0,len(r_array)):
        file.write('%.15e %.15e\n'%(r_array[i],rho_array[i])) # col1:r (AU), col2: surf. density (g/cm2)
    print('File generated!')
else: 
    print('Sorry, density profile generated does not reproduce total dust mass.')

sys.exit()




############################################################
# Case 0. Fictitious profile (gaussian gap)
def fictitious_density_profile_gaussian(r,sigma,mu):
    delta=1e-20 # Depletion factor (dimensionless)
    exp_value=-0.5*(r-mu)**2/sigma**2
    value_o=-0.99999999999999999999*(np.exp(exp_value)/(sigma*(2*np.pi)**0.5))+1
    sigma_0=1

    try:
        """
        if R_in<=r<d_min or d_max<r<=R_out:
            value=sigma_0*R_exp/r*np.exp(-r/R_exp)
            return value
        elif d_min<r<d_max :
            value=(sigma_0*R_exp/r*np.exp(-r/R_exp))*delta
            return value
        """
        if R_in<=r<=a_b:
            value=sigma_0*R_exp/r*np.exp(-r/R_exp)
            return value
        #elif a_b<r<=29.0:
        #    value=1e-10
        #    return value
        elif a_b<r<=R_out:
            value=sigma_0*R_exp/r*np.exp(-r/R_exp)
            return value*value_o*sigma*(2*np.pi)**0.5
        else:
            raise ValueError
    except ValueError:
        return('r is not between R_in and R_out')                                    


############################################################
# Case 1. No gaps
def surface_density_no_gaps(r):
    
    ############################################################
    #
    # Returns the unperturbed density (that with no gaps) 
    # in g/cm2. The input is the radial distance in AU.
    #
    ############################################################

    exp_value=np.exp(-R_in/R_exp)-np.exp(-R_out/R_exp)
    sigma_0=M_dust_disk/(2.0*np.pi*R_exp**2*exp_value) * (cte.M_sun.value*1000/(cte.au.value*100)**2)
    value=sigma_0*R_exp/r*np.exp(-r/R_exp)
    return value


############################################################
# Case 2. One gap
def surface_density_one_gap(r,m,a):

    ############################################################
    #
    # Find the density surface of a disk with ONE planet inside.
    # The planet has mass 'm' (in solar masses) and it is located
    # at a distance 'a' (in AU) from the central star.
    # Output is the surface density (in g/cm^2) in the desired 
    # array of points 'r' (in AU).
    #
    ############################################################

    planet_b=Companion(m,a)
    K=(m/m_a)**2*(h_scale(m_a,a)/a)**-5*alpha**-1
    d_min_b=a-0.5*planet_b.gap_width()
    d_max_b=a+0.5*planet_b.gap_width()
    exp_value= np.exp(-R_in/R_exp)-np.exp(-R_out/R_exp)+\
               np.exp(-d_max_b/R_exp)-np.exp(-d_min_b/R_exp)+\
               (np.exp(-d_min_b/R_exp)-np.exp(d_max_b/R_exp))/(1+0.04*K)
    sigma_0=M_dust_disk/(2.0*np.pi*R_exp**2*exp_value) * (cte.M_sun.value*1000/(cte.au.value*100)**2)
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
    sigma_0=M_dust_disk/(2.0*np.pi*R_exp**2*exp_value) * (cte.M_sun.value*1000/(cte.au.value*100)**2)
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
