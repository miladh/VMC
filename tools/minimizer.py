# -*- coding: utf-8 -*-
"""
Created on Tue May 28 20:21:40 2013

@author: Milad H. Mobarhan
"""

from sys import argv
from pylab import *
import matplotlib.gridspec as gridspec

close("all")
filename = argv[1]

data = loadtxt(filename+ "/results", skiprows=1)
alpha = array(data[:,0])
beta = array(data[:,1])
energy  = array(data[:,2])
energyDer  = array(data[:,7])
steps = range(0,len(alpha))

#energy[117]=-14.48
#energy[76]=-14.48
#
#energyDer[76]=3.11789000e-03
#energyDer[117]=3.11789000e-03

alphaT = cumsum(alpha)/cumsum(ones(alpha.shape))
betaT = cumsum(beta)/cumsum(ones(beta.shape))
energyDerT = cumsum(energyDer)/cumsum(ones(energyDer.shape))
energyT = cumsum(energy)/cumsum(ones(energy.shape))

plt.figure(figsize=(15, 8))
plt.rc('font', size=11) 
G = gridspec.GridSpec(3, 2)


ax1=plt.subplot(G[0, 0])
plot(steps,alpha,label = "Exact")
plot(alphaT,'r',label = "Trailing Average")
ylabel(r'$\alpha$', fontsize=16)
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=2, mode="expand", borderaxespad=0.)
setp( ax1.get_xticklabels(), visible=False)


ax2=plt.subplot(G[1,0],sharex=ax1)
plot(steps,beta,label = "Exact")
plot(betaT,'r',label = "Trailing Average")
ylabel(r'$\beta$', fontsize=16)
setp( ax2.get_xticklabels(), visible=False)

ax3=plt.subplot(G[2, 0],sharex=ax1)
plot(steps,energy)
plot(energyT,'r')
ylabel(r'$E$', fontsize=16)
xlabel('Step')
#xlim(0,490)

#ax4=plt.subplot(G[3, 0],sharex=ax1)
#plot(steps,energyDer)
#plot(energyDerT,'r')
#ylabel(r'$\partial E/\partial \alpha$', fontsize=16)
#xlabel('Step')
#plt.tight_layout()


plt.subplot(G[:, 1])
plot(alphaT,beta,'-*')
xlabel(r'$\alpha$')
ylabel(r'$\beta$', fontsize=16)
#plt.tight_layout()


#plt.savefig(filename+"/minBe.pdf",bbox_inches='tight')





#figure()
#hist(data[:,7], bins=50)