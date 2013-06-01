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
import matplotlib.gridspec as gridspec


close("all")
defualtPath = "/home/milad/Dropbox/fys4411/code/vmc/DATA/onebodyDensity"

dataPaths=[]
for i in range(1,len(sys.argv)):
    path = join(defualtPath + argv[i])
    dataPaths.append(path)

i=0
Data=[]


#Plot 2d distribution
fig = plt.figure(figsize=(15, 8))
G = gridspec.GridSpec(3, 2)
#plt.subplot(1,2,1)
#xlabel(r'$r$',fontsize=20)
#ylabel(r'$r^2\rho(r)$',fontsize=20)
#xlim(0,5)
#labels = ("Pure", 'With Jastrow', "No Jastrow")


labels = ("Symmetric", 'Antisymmetric')


for paths in dataPaths:
    fileList=glob.glob1(paths,"*.bin")

    rawfilePaths =[]
    for files in fileList:
        rawfilePaths.append(join(paths,files))
    
    
    data = []
    for files in rawfilePaths:
        filename = join(dataPaths,files)
        data.append(fromfile(files,dtype=dtype(("float",3))))
    
    data = array(data)
    data = data.reshape(len(fileList)*len(data[0]),3)
    
    Data.append(data)
    x=Data[i][:,0]
    y=Data[i][:,1]
    z=Data[i][:,2]
    
    
#    r=zeros((len(x),1))
#    r=sqrt(x**2+y**2+z**2)
#    myBins=hist(r, bins=500, histtype='step',normed=True,label =labels[i]) 
    plt.subplot(G[:, 0])
    plt.xlim(-7,7)
    xlabel(r'$x$',fontsize=20)
    ylabel(r'$\rho(x)$',fontsize=20)
    myBins=hist(x, bins=500, histtype='step',normed=True,label =labels[i])
    
#    plt.subplot(G[1, 0])   
#    plt.xlim(-3,3)
#    xlabel(r'$y$',fontsize=20)
#    ylabel(r'$\rho(y)$',fontsize=20)
#    myBins=hist(y, bins=500, histtype='step',normed=True,label =labels[i])
#   
#    plt.subplot(G[2, 0])
#    plt.xlim(-3,3)
#    xlabel(r'$z$',fontsize=20)
#    ylabel(r'$\rho(z)$',fontsize=20)
#    myBins=hist(z, bins=500, histtype='step',normed=True,label =labels[i])
    i = i+1   
    
    print "done reading from file!"

#legend()


#Plot countour 2d 
Hbins=300
H, xedges, yedges = np.histogram2d(y, x, bins=(Hbins, Hbins),normed=True)
extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
H = ndimage.gaussian_filter(H, sigma=3.0, order= 0)

plt.subplot(G[:,1])
cont = plt.imshow(H, extent=extent)
plt.xlabel("$X$",fontsize=20)
plt.ylabel("$Y$",fontsize=20)
plt.axis([-6,6,-6 ,6])
plt.tight_layout()

cbar_ax = fig.add_axes([0.54, 0.032, 0.45, 0.02])
fig.colorbar(cont, cax=cbar_ax,orientation= "horizontal")




#figure size & labelsize
#fig = plt.figure()
#plt.rc('font', size=11)

#Plot 3d distribution
#xmesh = linspace(xedges[0],xedges[-1],Hbins)
#ymesh = linspace(yedges[0],yedges[-1],Hbins)
#Y, X = np.meshgrid(ymesh,xmesh)
#
#ax = fig.add_subplot(2, 2, 1,projection='3d')
#surf = ax.plot_surface(X, Y, H, rstride=1, cstride=1,
#                       cmap=cm.Blues,linewidth=0.05)
#ax.set_ylabel(r"$X$",fontsize=20)
#ax.set_xlabel(r"$Y$",fontsize=20)
#ax.set_zlabel(r"$\rho$",fontsize=20)
#
#
#fig.tight_layout()
#ax.view_init(40,50)
#maxValueCoorIndecies = where(H == H.max())
#cset = ax.contour(X, Y, H, zdir='x', offset=xedges[0]*1.05,levels=[0])
#cset = ax.contour(X, Y, H, zdir='y', offset=yedges[0]*1.05,levels=[1.0])
#
#
#
##colorbar
##fig.colorbar(surf, shrink=0.5, aspect=5)

