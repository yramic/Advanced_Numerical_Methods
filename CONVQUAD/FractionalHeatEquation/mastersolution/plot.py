import pandas as pd
import matplotlib.pyplot as plt
from sys import argv

input = str(argv[1])
output= str(argv[2])
# Load the data into a DataFrame
columns = ["M","time_MOT","time_TOEP","time_ASAO"]
df = pd.read_csv(input, usecols=columns)

# Create a single plot for all values of 'c'
plt.figure()
plt.title('Runtimes')

# Iterate through each group and plot the data
plt.plot(df.M, df.time_MOT, 'o-', label="MOT")
plt.plot(df.M, df.time_TOEP, 'o-', label="TOEP")
plt.plot(df.M, df.time_ASAO, 'o-', label="ASAO")
plt.plot(df.M, df.M/1000, '-.', color="black" ,label=r"$(\mathcal{O}(M))$")


plt.xlabel('M')
plt.ylabel("time[s])")
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig(output)
