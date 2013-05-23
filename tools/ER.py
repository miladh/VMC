# -*- coding: utf-8 -*-
"""
Created on Sun May 19 01:16:01 2013

author: Milad H. Mobarhan
"""

from sys import argv
from pylab import *

filename = argv[1]
   
data = loadtxt(filename, skiprows=2)
blokSize = data[0,:]
sigma = data[2,:]

plot(blokSize, sigma, "*")
grid()
xlabel('Block size')
ylabel('$\sigma$', fontsize=16)
show()