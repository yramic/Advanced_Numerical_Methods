
from numpy import *
from pylab import *
import sys

errors = loadtxt("FSG_errors.txt")
numberOfElements = loadtxt("FSG_N.txt")
print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)

loglog(numberOfElements, errors, "-o")

xlabel("discretization parameter N")
ylabel("$||rho_{N}-rho||$")
grid("on")
show()

print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)
