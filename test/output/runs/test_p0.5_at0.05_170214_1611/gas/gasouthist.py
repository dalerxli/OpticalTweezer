import math
import numpy as np

data = np.genfromtxt("out_run_2200",usecols=3)

histogram,bins = np.histogram(data,bins=100,normed=True)
vsqdsum = np.mean(data)
# for i in range(0,len(data)):
    # vsqdsum += data[i]

vsqdsum = vsqdsum*0.1/(3.)
print(vsqdsum)
outfile = open("outtemphist.dat","w")

for i in range(0,len(histogram)):
    outfile.write(str(bins[i]) + "\t" + str(histogram[i]) + "\n")
