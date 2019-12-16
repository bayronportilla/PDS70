import numpy as np
from astropy import constants as cte
from astropy import units
from astropy.io import fits
import matplotlib.pyplot as plt
import sys

hdulist=fits.open("PDS_70_2016-03-26_QPHI_CROPPED.fits")
data=hdulist[0].data
#hdulist=fits.open("PDS70_cont-final.fits")
#data=hdulist[0].data[0][0]

class Pixel:
    # constructor
    def __init__(self,x,y,z):
        self.x=x # atribute
        self.y=y # atribute
        self.z=z # atribute

    def gotocenter(self,xmax,ymax): # method
        self.x=int(self.x-0.5*xmax)
        self.y=int(self.y-0.5*ymax)
        pass
        
    def deproject(self,pa,inc): # method
        pa=(pa*units.deg).to(units.rad).value
        inc=(inc*units.deg).to(units.rad).value
        # The following matrix rotate a vector by an angle 
        # theta, counterclockwise if |angle|>0
        M=np.array([[np.cos(pa),-np.sin(pa)],[np.sin(pa),np.cos(pa)]])
        xxx=self.x
        yyy=self.y
        xx=np.dot(M,np.array([xxx,yyy]))[0]
        yy=np.dot(M,np.array([xxx,yyy]))[1]
        x=xx*np.cos(inc)
        y=yy
        self.x=int(round(x))
        self.y=int(round(y))
        pass

    def gotocorner(self,xmin,ymin): # method
        self.x=int(self.x+abs(xmin))
        self.y=int(self.y+abs(ymin))
        
    #def originatcorner(self,xmin):

parray=[]
for i in range(0,data.shape[0]):
    for j in range(0,data.shape[1]):
        p=Pixel(i,j,data[i][j])
        parray.append(p)

xmax=max((p.x) for p in parray)
ymax=max((p.y) for p in parray)

for p in parray:
    p.gotocenter(xmax,ymax)
    p.deproject(158.6,49.7)

xmin=min((p.x) for p in parray)
ymin=min((p.y) for p in parray)

for p in parray:
    p.gotocorner(xmin,ymin)

xmax=max((p.x) for p in parray)
ymax=max((p.y) for p in parray)

data_deprojected=np.zeros((xmax+1,ymax+1))

for p in parray:
    data_deprojected[p.x][p.y]=p.z

hdu=fits.PrimaryHDU(data_deprojected)
#hdu1=fits.HDUlist([hdu])
hdu.writeto("Jband_deprojected.fits")
sys.exit()

plt.imshow(data_deprojected,origin="lower")
#plt.imshow(data_deprojected)
plt.axes().set_aspect("equal")
plt.show()



    
    




    





    





        



