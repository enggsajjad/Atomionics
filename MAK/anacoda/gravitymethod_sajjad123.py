#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# ## Point Mass<br>
# we derived the vertical component of the gravitational acceleration due to a point mass with mass m as:
# 
# $\Delta_g = \frac {Gmz}{(x^2+z^2)^{3/2}}$
# 

# In[2]:


#Assuming a ponit mass 
start = 500
N = 4000.0
G = 6.76e-11
x = np.linspace(-(start),(start),(start*2+1))
m = 1550.0
Zi = 175             
di = -1             

delta_g=0.0*x
for tau in x:
    delta_g[x==tau] = ((G*m*Zi / (((tau-di)**2+Zi**2)**(3/2))))*10**(6)

plt.plot(x,delta_g)
plt.xlabel('distance (x)',fontsize=15)
plt.ylabel('Gravitational Accel.$\Delta_g$',fontsize=15)
plt.title('Ponit Mass Gravitational Acceleration')


# ## Complex Body with Multiple Point Masses<br>
# We can approximate the body with complex shape as a distribution of point masses. 
# The gravitational attraction of the body is then nothing more than the sum of the gravitational attractions of all of 
# the individual point masses as illustrated below.
# If there are N point masses, this equation can be written more compactly in mathematical notation, this sum can be written as:
# 
# $\Delta_g = \sum_{i=0}^N \frac {Gmz_i}{((x-d_i)^2+z_i^2)^{3/2}}$
# 
# where z represents the depth of burial of each point mass, d represents the horizontal position of each point mass, 
# and x represents the horizontal position of the observation point. Only the first three terms have been written in 
# this equation. There is, in actuality, one term in this expression for each point mass.
# 
# ![title](image.jpg)

# In[3]:


#Assuming a mass of 3x3point square
start = 200
N = 1000.0

G = 6.76e-11

x = np.linspace(-(start),(start),(start*2+1))
m = 1.0
#assum these values for Zi
#Zi = np.array([4.0,4.0,4.0,6.0,6.0,6.0,8.0,8.0,8.0])
#Zi = np.array([4.0,8.0,14.0,16.0,9.0,10.0,18.0,28.0,15.0])
Zi = np.random.rand(9,1)*10+10
#assum these values for di
#di = np.array([-1.0,0.0,1.0,-1.0,0.0,1.0,-1.0,0.0,1.0])
#di = np.array([-5.0,0.0,1.0,-10.0,4.0,21.0,-31.0,4.0,51.0])
di = np.random.rand(9,1)*100-50
#%%
delta_g=0.0*x
delta_g1=0.0*x
delta_g2=0.0*x
delta_g3=0.0*x
delta_g4=0.0*x
delta_g5=0.0*x
delta_g6=0.0*x
delta_g7=0.0*x
delta_g8=0.0*x
delta_g0=0.0*x

for tau in x:
    delta_g[x==tau] = np.sum(G*m*Zi / (((tau-di)**2+Zi**2)**(3/2)))
    delta_g1[x==tau] = (G*m*Zi[1] / (((tau-di[1])**2+Zi[1]**2)**(3/2)))
    delta_g2[x==tau] = (G*m*Zi[2] / (((tau-di[2])**2+Zi[2]**2)**(3/2)))
    delta_g3[x==tau] = (G*m*Zi[3] / (((tau-di[3])**2+Zi[3]**2)**(3/2)))
    delta_g4[x==tau] = (G*m*Zi[4] / (((tau-di[4])**2+Zi[4]**2)**(3/2)))
    delta_g5[x==tau] = (G*m*Zi[5] / (((tau-di[5])**2+Zi[5]**2)**(3/2)))
    delta_g6[x==tau] = (G*m*Zi[6] / (((tau-di[6])**2+Zi[6]**2)**(3/2)))
    delta_g7[x==tau] = (G*m*Zi[7] / (((tau-di[7])**2+Zi[7]**2)**(3/2)))
    delta_g8[x==tau] = (G*m*Zi[8] / (((tau-di[8])**2+Zi[8]**2)**(3/2)))
    delta_g0[x==tau] = (G*m*Zi[0] / (((tau-di[0])**2+Zi[0]**2)**(3/2)))
    
plt.plot(x,delta_g0)
plt.plot(x,delta_g1)
plt.plot(x,delta_g2)
plt.plot(x,delta_g3)
plt.plot(x,delta_g4)

plt.plot(x,delta_g5)
plt.plot(x,delta_g6)
plt.plot(x,delta_g7)
plt.plot(x,delta_g8)

plt.plot(x,delta_g,'r--')
plt.xlabel('distance (x)',fontsize=15)
plt.ylabel('Gravitational Accel.$\Delta_g$',fontsize=15)
plt.title('Multiple Ponit Mass Gravitational Acceleration')


# In[4]:
#%%

plt.plot(di,Zi,'r*')
plt.xlim([-start, start])
plt.ylim([0, 5+np.ceil(max(Zi))])
plt.xlabel('distance (d)',fontsize=15)
plt.ylabel('Depth (Z)',fontsize=15)
plt.title('Rock Multiple Masses')

