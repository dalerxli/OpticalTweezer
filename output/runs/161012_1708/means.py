import numpy as np

data_in = np.loadtxt("gasTempData.dat",usecols=(0,))
data_out = np.loadtxt("gasTempData.dat",usecols=(1,))

print("In: " + str(np.mean(data_in)) + " +- " + str(np.std(data_in)))
print("Out: " + str(np.mean(data_out)) + " +- " + str(np.std(data_out)))
