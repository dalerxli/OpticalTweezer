import math
import numpy as np

data = np.genfromtxt("out_run_15400",skiprows=100,usecols=3)

histogram,bins = np.histogram(data,bins=100,normed=True)

outfile = open("outtemphist.dat","w")

for i in range(0,len(histogram)):
    outfile.write(str(bins[i]) + "\t" + str(histogram[i]) + "\n")
