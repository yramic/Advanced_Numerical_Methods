
from numpy import *
from pylab import *
import sys

errors = loadtxt("integrateTiSing_errors.txt")
numberOfElements = loadtxt("integrateTiSing_N.txt")

loglog(numberOfElements, errors, "-o")

xlabel("quadrature points")
ylabel("error")
grid("on")
show()
