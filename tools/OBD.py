# -*- coding: utf-8 -*-
"""
Created on Sat Feb  9 15:38:09 2013

@author: Milad H. Mobarhan 
"""

from mpl_toolkits.mplot3d import Axes3D
from pylab import*
from sys import argv
from os.path import join
import glob
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage


close("all")
#dataPath = argv[1]
#fileList=glob.glob1(dataPath,"*.mat")
#
#data=[]
#for files in fileList:
#    filename= join(dataPath,files)
#    data.append(loadtxt(filename))
#
#data = array(data)
#data = data.reshape(len(fileList)*len(data[0]),3)
#
#x=data[:,0]
#y=data[:,1]
#z=data[:,2]
#r=zeros((len(x),1))
#
#for i in range(len(x)):
#    r[i]=sqrt(x[i]**2+y[i]**2+z[i]**2)


#Plot 2d distribution
myBins=hist(r, bins=100, histtype='step',normed=True) 
xlabel(r'$r$',fontsize=20)
ylabel(r'$\rho(r)$',fontsize=20)
grid()




#Plot countour 2d 
plt.figure()
H, xedges, yedges = np.histogram2d(y, x, bins=(100, 100),normed=True)
extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
#H = ndimage.gaussian_filter(H, sigma=3.0, order= 0)


plt.imshow(H, extent=extent)
plt.xlabel("$X$",fontsize=20)
plt.ylabel("$Y$",fontsize=20)
plt.colorbar()
plt.show()


#figure size & labelsize
fig = plt.figure()
plt.rc('font', size=11) 

#Plot 3d distribution
Y, X = np.meshgrid(linspace(yedges[0],yedges[-1],100),linspace(xedges[0],xedges[-1],100))

ax = fig.add_subplot(1, 1, 1,projection='3d')
surf = ax.plot_surface(X, Y, H, rstride=1, cstride=1,cmap=cm.Blues,linewidth=0)
ax.set_ylabel(r"$X$",fontsize=20)
ax.set_xlabel(r"$Y$",fontsize=20)
ax.set_zlabel(r"$\rho$",fontsize=20)

fig.tight_layout()
ax.view_init(40,50)
maxValueCoorIndecies = where(H == H.max())
cset = ax.contour(X, Y, H, zdir='x', offset=xedges[0]*1.05,levels=[0])
cset = ax.contour(X, Y, H, zdir='y', offset=yedges[0]*1.05,levels=[2.3])



#colorbar
fig.colorbar(surf, shrink=0.5, aspect=5)

