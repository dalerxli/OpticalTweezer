import numpy as np
import math

temps = np.genfromtxt("gasTempData.dat",usecols=(0,1),skiprows=22)
# for i in range(0,len(temps)):
    # if math.isnan(temps[i][1]):
        # print(i)
        # # del temps[i][1]
# print(temps[:,0])
# print(temps[:,1])
print("Mean inTemp: " + str(np.mean(temps[:,0])))
print("Mean outTemp: " + str(np.mean(temps[:,1])))
print("Mean difference " + str(np.mean(temps[:,1])-np.mean(temps[:,0])))

for i in range(0,len(temps)):
    print(temps[i][1] - temps[i][0])




