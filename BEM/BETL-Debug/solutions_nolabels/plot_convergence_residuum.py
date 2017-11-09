
from numpy import *
from pylab import *
import sys

numberOfElements = loadtxt("BETL-Debug_levels.txt")
rD = loadtxt("BETL-Debug_rDnorm.txt")
rN = loadtxt("BETL-Debug_rNnorm.txt")

line1, = loglog(numberOfElements, rD, "-o")
line2, = loglog(numberOfElements, rN, "-*")

plt.legend([line1, line2], [r'$\rho_D$', r'$\rho_N$'])
xlabel("Number of elements")
ylabel(r"$||\rho||_{\infty}$")
grid("on")
show()

print polyfit(log(numberOfElements), log(rD), 1)

print polyfit(log(numberOfElements), log(rN), 1)
