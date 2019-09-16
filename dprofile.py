import numpy as np
import matplotlib.pyplot as plt
import sys

xval=np.linspace(-1,1,1000)
mu=0.0
sigma=0.3
density=[]

# Gaussian depletion factor
def gdepf(x_array,mu,sigma):
    value=1/(sigma*(2*np.pi)**0.5) * np.exp(-0.5*((x_array-mu)/sigma)**2)
    return value

def delta(x_array,mu,sigma):
    ymax=max(gdepf(x_array,mu,sigma))
    return gdepf(x_array,mu,sigma)/ymax



def dbefore(x_array):
    value=x_array/x_array
    return value



def dafter(x_array,mu,sigma):
    da=[]
    for i in range(0,len(x_array)):
        if x_array[i]<=0.0:
            da.append(dbefore(x_array[i]))
        else:
            da.append(dbefore(x_array[i])*delta(x_array,mu,sigma)[i])
    return da


    
#density=np.array(density)
    

fig,((ax_1),(ax_2),(ax_3))=plt.subplots(1,3,figsize=(14,8))
ax_1.plot(xval,delta(xval,mu,sigma),".")
ax_2.plot(xval,dbefore(xval),".")
ax_3.plot(xval,dafter(xval,mu,sigma),".")
ax_1.set_xlim(-1,1)
ax_2.set_xlim(-1,1)
ax_3.set_xlim(-1,1)
plt.show()
