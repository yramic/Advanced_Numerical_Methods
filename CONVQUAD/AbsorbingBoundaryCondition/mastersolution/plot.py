import numpy as np
import matplotlib.pyplot as plt
from sys import argv

# file names
data = str(argv[1])
output = str(argv[2])

# read solution data
d = np.loadtxt(data)

# prepare grid
x = np.linspace(0, 1, d.shape[0])
t = np.linspace(0, 1, d.shape[1])
X, T = np.meshgrid(x, t)
Z = np.transpose(d)

fig, ax = plt.subplots()
CS = ax.contour(X, T, Z)
# for visual clarity
manual_locations = [
    (0.4, 0.05), (0.5, 0.1), (0.6, 0.2), (0.5, 0.3), (0.2, 0.4), (0.4, 0.6), (0.5, 0.7), (0.6, 0.75), (0.7, 0.85), (0.8, 0.9), (0.9, 0.95)]
ax.clabel(CS, inline=True, fontsize=10, manual=manual_locations)
ax.set_title('$u(x,t)$')
ax.set_xlabel("$x$")
ax.set_ylabel("$t$")

plt.savefig(output, dpi=150)
