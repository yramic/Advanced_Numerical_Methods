
from numpy import *
from pylab import *
import sys

errors = loadtxt("AR_errors.txt")
numberOfElements = loadtxt("AR_N.txt")
print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)

loglog(numberOfElements, errors, "-o")

xlabel("discretization parameter N")
ylabel("$||u_{N}-u||$")
grid("on")
show()

print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)
