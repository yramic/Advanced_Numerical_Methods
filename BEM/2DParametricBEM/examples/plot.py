import numpy as np
import matplotlib.pyplot as plt

fname = "log_out.txt"
data = np.loadtxt(fname)

x = data[:,0]
y = data[:,1]
order = data[:,3]

fig,ax1 = plt.subplots()
ax1.plot(x,y,'b')

ax1.set_xlabel('y')
ax1.set_ylabel('$\int_{-1}^{1} \log{(x-y)} dx$')
#ax1.legend(['f(y)'])
ax1.tick_params('y',colors='b')
ax2 = ax1.twinx()
ax2.plot(x,order,'r')
ax2.set_ylabel('Gauss quadrature order')
#ax2.legend(['Order'])
ax2.tick_params('y',colors='r')
plt.title("Plot of $\int_{-1}^{1} \log{(x-y)} dx$ vs y")
plt.show()
