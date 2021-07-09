# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 11:04:17 2021

@author: Dell
"""

import numpy as np
import matplotlib.pyplot as plt 
from scipy.stats import gaussian_kde 
from mpl_toolkits import mplot3d

rng = np.arange(-1000,1025,25, dtype = 'f8')

#part 1
y =  8*(rng**3) + 12*rng
fig, ax = plt.subplots()
ax.scatter(rng,y, s = 5, c = 'red', label = '8x^3 + 12x (scatter)')
ax.plot(rng,y, label = '8x^3 + 12x (line)')
ax.set_ylim(-10000000000,10000000000)
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.set_title('x vs f(x)')
ax.legend()
plt.savefig("hw1a.png")
plt.show()

#part 2 8*(n-10*(m-25))3 – 12*(n-10*(m-25))
ax2 = plt.axes()
rng2 = np.arange(0,51)
color = [1,0,0,1]
for i in rng2:
    y2 = 8*((i-25)*(rng-10))**3 - 12*((i-25)*(rng-10))
    if i <= 10:
        color[1] += 10/111
    elif i<=20:
        color[0] -= 10/111
    elif i<=30:
        color[2]+=10/111
    elif i<=40:
        color[1] -=10/111
    elif i<=50:
        color[0] += 10/111
    tup = (color[0],color[1],color[2],color[3])
    ax2.plot(rng,y2,c = tup)
ax2.set_xlabel('x')
ax2.set_ylabel('f(x)')
ax2.set_title('x vs f(x)')
plt.savefig("hw1b.png")
plt.show()

#part 3
#fig2, ax3 = plt.subplots(figsize = (10,4))
fig3, ax4 = plt.subplots(figsize = (10,4))
filenames = ['data1.dat','data2.dat','data3.dat','data4.dat']
cmaps = ['spring', 'summer', 'autumn', 'winter']
for count,file in enumerate(filenames):
    arr = np.genfromtxt(file,skip_header = 1, usecols = (8,12), dtype = float, max_rows = 2000)
    x = arr[:, 0]#RGYR_FOLDED
    y = arr[:, 1]#THET
    xy = np.vstack((x,y))
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x,y,z = x[idx], y[idx], z[idx] 
#    x1 = ax3.scatter(x,y, c = z, cmap = cmaps[count])
#    fig2.colorbar(x1, shrink = 0.6)
    ax4=plt.axes(projection='3d')
    ax4.scatter3D(x, y, z, c = z, cmap=cmaps[count])
#ax3.xlabel('Rg ($\AA$)', fontsize=20)
#ax3.ylabel('THET ', fontsize=20)
#ax3.set_title('probability density')



#    ax.set_xlabel('RMSD ($\AA$)', fontsize=16)
#    ax.set_ylabel('Rg ($\AA$)', fontsize=16)
#    ax.set_zlabel('P(RMSD,Rg) ($\AA^{-2}$)', fontsize=16)



