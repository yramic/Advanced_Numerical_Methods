
from numpy import *
from pylab import *
import sys

errors = loadtxt("GLQR_errors.txt")
errors2 = loadtxt("WLQR_errors.txt")
n = loadtxt("QR_N.txt")

line1, = loglog(n, errors, "-o", color="blue")
line2, = loglog(n, errors2, "-*", color="red")
plt.legend([line1, line2], ['Gauss-Laguerre', 'Log-Weighted'])

xlabel("discretization parameter N")
ylabel("$error$")
grid("on")
show()

print diff(log(errors)) / diff (log(numberOfElements))
print polyfit(log(numberOfElements), log(errors), 1)
