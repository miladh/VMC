from sys import argv
from pylab import *
close("all")
filename = argv[1]
   
data = loadtxt(filename+"/blocking.mat", skiprows=2)
blokSize = data[0,:]
sigma = data[2,:]


plot(blokSize, sigma,'*')
xlabel('Block size')
ylabel('$\sigma(m_n)$', fontsize=16)
plt.tight_layout()





#x=array(0.00769,0.091,0.057,0.131,0.164,0.146,0.15)
#c=logspace(1,7,7)