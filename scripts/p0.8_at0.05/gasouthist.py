import math
import numpy as np

data = np.genfromtxt("p0.8_at0.05_170207_0344_gasTempData.dat",skiprows=100,usecols=1)

histogram,bins = np.histogram(data,bins=100,normed=True)

outfile = open("outtemphist.dat","w")

for i in range(0,len(histogram)):
    outfile.write(str(bins[i]) + "\t" + str(histogram[i]) + "\n")
