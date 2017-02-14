import math
import numpy as np

data="plotfile_clean.dat"

temdata = np.genfromtxt(data,usecols=7)
timpdata = np.genfromtxt(data,usecols=5)
tsurdata = np.genfromtxt(data,usecols=3)

alpha = []

for i in range(0,len(temdata)):
    temp = (temdata[i]-timpdata[i])/(tsurdata[i]-timpdata[i])
    alpha.append(temp)
    print(temp)

print("")
print(np.mean(alpha))
