# -*- coding: utf-8 -*-
"""
Created on Sun May 19 01:16:01 2013

author: Milad H. Mobarhan
"""
from sys import argv
from pylab import *

close("all")


#defualtPath = "/home/milad/Dropbox/fys4411/code/vmc/DATA/STminimizing/"
defualtPath = "/home/milad/Dropbox/fys4411/code/vmc/"

filename = argv[1]   
data = loadtxt(defualtPath+filename+ "/results", skiprows=1)
R = data[:,6]
energy = data[:,2]
alpha =data[:,0]
beta = data[:,1]


#filename = argv[2]   
#data = loadtxt(defualtPath+filename+ "/results", skiprows=1)
#R1 = data[:,6]
#energy1 = data[:,2]
#alpha1 =data[:,0]
#beta1 = data[:,1]

energyT = cumsum(energy)/cumsum(ones(energy.shape))


fig = plt.figure(figsize=(10, 8))
plt.rc('font', size=15) 
plot(R, energy, "--",linewidth=3,label=r"$\Psi_+$")
plot(R, energyT)
#plot(R1, energy1, "--",linewidth=3,label=r"$\Psi_-$")
xlabel('$R$',fontsize=20)
ylabel(r'$\langle E_L \rangle$', fontsize=20)
#xlim(0.5,3)
fig.tight_layout()
legend()






#fig = plt.figure(figsize=(15, 10))
#ax1=plt.subplot(2,1,1)
#plt.rc('font', size=18) 
#plot(R, alpha, "--",linewidth=3, label=r"$\Psi_+$")
#plot(R1, alpha1, "--",linewidth=3, label=r"$\Psi_-$")
#ylabel(r'$\alpha$', fontsize=20)
#setp( ax1.get_xticklabels(), visible=False)
#xlim(0.5,3)
#legend()
#
#ax2=plt.subplot(2,1,2)
#plot(R, beta, "--",linewidth=3, label=r"$\Psi_+$")
#plot(R1, beta1, "--",linewidth=3,label=r"$\Psi_--$")
#xlabel('$R$',fontsize=20)
#ylabel(r'$\beta$', fontsize=20)
#xlim(0.5,3)
#fig.tight_layout()
