import numpy as np
import matplotlib.pyplot as plt

fname1 = "./../build/sl_circle_test_eigs.txt"
data = np.loadtxt(fname1)
# Getting the unique eigenvalues
eigs_unique = data[1::2]
eigs_unique = np.flip(eigs_unique,0)
test = 1/eigs_unique
plt.figure()
plt.plot(test)
plt.ylabel("$eigenvalue^{-1}$")
plt.show()
