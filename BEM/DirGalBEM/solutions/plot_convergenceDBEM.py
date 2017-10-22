
from numpy import *
from pylab import *
import sys

errors = loadtxt("DBEM1stK_errors.txt")
errors2 = loadtxt("DBEM2ndK_errors.txt")
numberOfElements = loadtxt("BEM_N.txt")
print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)

line1, = loglog(numberOfElements, errors, "-o", label='Line 1')
line2, = loglog(numberOfElements, errors2, "--", label='Line 2')
plt.legend([line1, line2], ['DBEM 1stKind', 'DBEM 2ndKind'])

xlabel("Number of dofs")
ylabel("$||u_{N}-u||_V$")
grid("on")
show()

print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)
