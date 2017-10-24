
from numpy import *
from pylab import *
import sys

errors = loadtxt("DBEM2ndK_errors.txt")
errors2 = loadtxt("DBEM2ndK_L2errors.txt")
numberOfElements = loadtxt("BEM_N.txt")
print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)

line1, = loglog(numberOfElements, errors, "-o", color='blue', label='Line 1')
line2, = loglog(numberOfElements, errors2, "--*", color='red', label='Line 2')
plt.legend([line1, line2], ['$||u_{N}-u||_V$', '$||u_{N}-u||_{L2}$'])
plt.title('DBEm Second-kind')

xlabel("Number of dofs")
ylabel("error")
grid("on")
show()

print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)
