import numpy as np
import matplotlib.pyplot as plt

fname1 = "sl_circle_test_eigs.txt"
data = np.loadtxt(fname1)
eigs_unique = data[1::2]
#test = test**2
eigs_unique = np.flip(eigs_unique,0)
test = 1/eigs_unique
#l = len(eigs_unique)
#x = np.linspace(1,l,l)
#print(x)
#print(1/eigs_unique)
plt.figure()
plt.plot(test)
plt.ylabel("$eigenvalue^{-1}$")
#plt.plot(1/x,'r')
plt.show()
