import numpy as np
L_sun=3.828e26
R_sun=6.9551e8 
L_b=1.6e-4*L_sun
L_c=1.6e-4*L_sun
R_b=1.26*R_sun
R_c=1.26*R_sun
m_b=7.5*(1.898e27/1.989e30) 
m_c=8.0*(1.898e27/1.989e30) 
sigma=5.67e-8

############################################################
# Source temperature
T_b=(L_b/(4*np.pi*sigma*R_b**2))**(1/4.)
T_c=(L_c/(4*np.pi*sigma*R_c**2))**(1/4.)

print('T_b: ',T_b)
print('T_c: ',T_c)
print('L_b: ',L_b/L_sun)
print('L_c: ',L_c/L_sun)
print('m_b: ',m_b)
print('m_c: ',m_c)

