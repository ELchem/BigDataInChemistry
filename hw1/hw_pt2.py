# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 19:17:05 2021

@author: Dell
"""
import numpy as np
import matplotlib.pyplot as plt 
from scipy.stats import gaussian_kde 
from mpl_toolkits import mplot3d

def get_all_figures():
    return [plt.figure(i) for i in plt.get_fignums()]

#part 3
fig2, ax3 = plt.subplots(figsize = (10,4))
fig = plt.figure(figsize = (7,4))
ax4 = fig.add_subplot(111, projection='3d')
plt.tight_layout()

filenames = ['data1.dat','data2.dat','data3.dat','data4.dat']
cmaps = ['spring', 'summer', 'autumn', 'winter']
x_all = np.empty(0)
y_all = np.empty(0)
xy_all = np.empty((0,0))

for count,file in enumerate(filenames):
    arr = np.genfromtxt(file,skip_header = 1, usecols = (8,12), dtype = float, max_rows = 2000)
    x = arr[:, 0]#RGYR_FOLDED
    y = arr[:, 1]#THET
    xy = np.vstack((x,y))
    x_all = np.append(x_all,x)
    y_all = np.append(y_all,y)
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x,y,z = x[idx], y[idx], z[idx] 
    x1 = ax3.scatter(x,y, c = z, cmap = cmaps[count])
    fig2.colorbar(x1, shrink = 0.6)
    ax4.scatter(x, y, z, c = z, cmap=cmaps[count])

#labels for 2d  
ax3.set_xlabel('Rg ($\AA$)', fontsize=20)
ax3.set_ylabel('THET ', fontsize=20)
ax3.set_title('probability density', fontsize = 20)
#labels for 3d
ax4.set_xlabel('Rg ($\AA$)', fontsize=20)
ax4.set_ylabel('THET ', fontsize=20)
ax4.set_zlabel('prob. dens.',fontsize = 20)
ax4.set_title('probability density', fontsize = 20)

#plot all datasets combined
xy_all = np.vstack((x_all,y_all))
fig3, ax5 = plt.subplots(figsize = (10,4))
z_all = gaussian_kde(xy_all)(xy_all)
idx = z.argsort()
x,y,z = x_all[idx], y_all[idx], z_all[idx] 
x2 = ax5.scatter(x_all, y_all, c = z_all, cmap = cmaps[count])

#labels for 2d combined scatter;
fig3.colorbar(x1, shrink = 0.6)
ax5.set_xlabel('Rg ($\AA$)', fontsize=20)
ax5.set_ylabel('THET ', fontsize=20)
ax5.set_title('probability density', fontsize = 20)

#save figs
plt.savefig('scatter2d_combined.png', dpi=300, bbox_inches='tight')
plt.figure(fig.number)
plt.savefig('scatter3d.png', dpi=300, bbox_inches='tight')
plt.figure(fig2.number)
plt.savefig('scatter2d.png', dpi=300, bbox_inches='tight')
