from mcmax3d_analysis.mcmax3d_deprojection import projected_cartesian_coordinates
from mcmax3d_analysis.mcmax3d_deprojection import deproject
import numpy as np
import matplotlib.pyplot as plt
import sys

pos_b=[96.8/1000,-147.9/1000]
pos_c=[-233.7/1000,28.8/1000]
phi=158.6
theta=49.7
d=113.4
print(deproject(pos_b[0],pos_b[1],theta,phi,d))
print(deproject(pos_c[0],pos_c[1],theta,phi,d))



sys.exit()
#print(projected_cartesian_coordinates(b=pos_b,c=pos_c,d=113.4))
x_array=np.linspace(-1,1,500)
x_array=list(x_array)
y_plus=[circle(i)[0] for i in x_array]
y_minus=[circle(i)[1] for i in x_array]
for i in range(0,len(x_array)):
    x_array.append(x_array[i])
y_array=y_minus+y_plus
x_array=np.array(x_array)
y_array=np.array(y_array)

xx_array=[steps(i,j,phi,theta)[0] for i,j in zip(x_array,y_array)]
yy_array=[steps(i,j,phi,theta)[1] for i,j in zip(x_array,y_array)]
xxx_array=[steps(i,j,phi,theta)[2] for i,j in zip(x_array,y_array)]
yyy_array=[steps(i,j,phi,theta)[3] for i,j in zip(x_array,y_array)]

x_d=[back(i,j,phi,theta)[0] for i,j in zip(xxx_array,yyy_array)]
y_d=[back(i,j,phi,theta)[1] for i,j in zip(xxx_array,yyy_array)]

fig=plt.figure()
ax=plt.axes()
ax.plot(x_array,y_array,'.')
ax.plot(xx_array,yy_array,'+')
ax.plot(xxx_array,yyy_array,'*')
ax.plot(x_d,y_d,'-')
ax.set_aspect('equal')
plt.show()
#print(y_array)



