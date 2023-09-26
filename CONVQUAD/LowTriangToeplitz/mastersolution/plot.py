import numpy as np
import matplotlib.pyplot as plt
from sys import argv

rt_mul = str(argv[1])
rt_lse = str(argv[2])
output = str(argv[3])
# read runtime data of matrix multiplication from text file
with open(rt_mul, "r") as file:
    d_m = [[float(x) for x in line.split(' ')] for line in file]
d_m = np.array(d_m)

# data for O(N^i), start and end adjusted for clearity
N1 = np.logspace(-6, np.log10(d_m[-1,0]/d_m[0,0])-6, d_m.shape[0])
N2 = np.logspace(-4.5, 2*np.log10(d_m[-1,0]/d_m[0,0])-4.5, d_m.shape[0])
N3 = np.logspace(-3, 3*np.log10(d_m[-1,0]/d_m[0,0])-3, d_m.shape[0])

fix, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,4))
ax1.loglog(d_m[:, 0], d_m[:, 1], 'b^-', markerfacecolor='none', label='matrix-matrix')
ax1.loglog(d_m[:, 0], d_m[:, 2], 'kx-', markerfacecolor='none', label='matrix-vector')
ax1.loglog(d_m[:, 0], d_m[:, 3], 'ro-', markerfacecolor='none', label='vector-vector')
ax1.loglog(d_m[:, 0], N3, 'b^:', markerfacecolor='none', label='$O(N^3)$')
ax1.loglog(d_m[:, 0], N2, 'kx:', markerfacecolor='none', label='$O(N^2)$')
ax1.loglog(d_m[:, 0], N1, 'ro:', markerfacecolor='none', label='$O(N^1)$')
ax1.legend()

# read runtime data of solving LSE from text file
with open(rt_lse, "r") as file:
    d_s = [[float(x) for x in line.split(' ')] for line in file]
d_s = np.array(d_s)

# data for O(N^i), start and end adjusted for clearity
N1 = np.logspace(-3.5, np.log10(d_s[-1,0]/d_s[0,0])-3.5, d_s.shape[0])
N2 = np.logspace(-6, 2*np.log10(d_s[-1,0]/d_s[0,0])-6, d_s.shape[0])

ax2.loglog(d_s[:, 0], d_s[:, 1], 'kx-', markerfacecolor='none', label='triSolve')
ax2.loglog(d_s[:, 0], d_s[:, 2], 'ro-', markerfacecolor='none', label='ltpSolve')
ax2.loglog(d_s[:, 0], N2, 'kx:', markerfacecolor='none', label='$O(N^2)$')
ax2.loglog(d_s[:, 0], N1, 'ro:', markerfacecolor='none', label='$O(N^1)$')
ax2.legend()
plt.savefig(output, dpi=150)
