# -*- coding: utf-8 -*-
"""
Created on Sat Feb  9 15:38:09 2013

@author: Milad H. Mobarhan 
"""
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from sys import argv
from os.path import join
import scipy.ndimage as ndimage
plt.close("all")


#Read in Data
dataPath = argv[1]
filename = join(dataPath,"results")
data = np.loadtxt(filename,skiprows=1)

#Create grid
x = np.unique(data[:,0])
y = np.unique(data[:,1])
Z = data[:,2].reshape(x.size, y.size)
Y, X = np.meshgrid(y, x)


#figure size & labelsize
fig = plt.figure(figsize=(15,6))
plt.rc('font', size=11) 

#Make surface plot
ZInt = ndimage.gaussian_filter(Z, sigma=3.0, order=0)
ax = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax.plot_surface(X, Y, ZInt, rstride=1, cstride=1, 
       cmap=cm.Blues,linewidth = 0.1)
      
ax.set_ylabel(r"$\beta$",fontsize=20)
ax.set_xlabel(r"$\alpha$",fontsize=20)
ax.set_zlabel(r"$E$",fontsize=20)
fig.tight_layout()
ax.view_init(16,51)




#Make contour plot
ZInt = ndimage.gaussian_filter(Z, sigma=1.9, order=0)
plt.rcParams['contour.negative_linestyle'] = 'solid'
ax = fig.add_subplot(1, 2, 2)
ax.contourf(Y,X, ZInt, 500, cmap=cm.Blues)
C = ax.contour(Y, X, ZInt, 11, colors='black', linewidth=.1)

ax.clabel(C, fontsize=15)
ax.set_xlabel(r"$\beta$",fontsize=20)
ax.set_ylabel(r"$\alpha$",fontsize=20)
fig.tight_layout()

#colorbar
cbar_ax = fig.add_axes([0.01, 0.042, 0.45, 0.02])
fig.colorbar(surf, cax=cbar_ax,orientation= "horizontal")