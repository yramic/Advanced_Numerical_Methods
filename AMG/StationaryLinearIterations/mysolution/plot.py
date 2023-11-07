import pandas as pd
import matplotlib.pyplot as plt
from sys import argv

input = str(argv[1])
output= str(argv[2])
# Load the data into a DataFrame
df = pd.read_csv(input)

# Group data by 'c'
grouped = df.groupby('c')

# Create a single plot for all values of 'c'
plt.figure()
plt.title('Convergence Plot')

# Iterate through each group and plot the data
for c, group in grouped:
    n_values = group['n']
    lambda_values = group['lambda(X)']

    plt.plot(n_values, lambda_values, 'o-', label=f'c={c}')

plt.xlabel('n')
plt.ylabel(r"$\lambda(X)$")
plt.xscale('log')
plt.legend()
plt.savefig(output)
