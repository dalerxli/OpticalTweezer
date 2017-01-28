import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from subprocess import call

# proc = call("ls *.dat",shell=True)
datetime = "170123_2033_"
gasTempDataIn = np.genfromtxt(datetime+"gasTempData.dat",usecols=0,skip_header=100)
gasTempDataOut = np.genfromtxt(datetime+"gasTempData.dat",usecols=1,skip_header=100)
vCOMData_x = np.genfromtxt(datetime+"vCOMData.dat",usecols=0,skip_header=100)
vCOMData_y = np.genfromtxt(datetime+"vCOMData.dat",usecols=1,skip_header=100)
vCOMData_z = np.genfromtxt(datetime+"vCOMData.dat",usecols=2,skip_header=100)

vSqd = []
for i in range(0,len(vCOMData_x)):
    vSqd.append(vCOMData_x[i]*vCOMData_x[i]+vCOMData_x[i]*vCOMData_x[i]+vCOMData_x[i]*vCOMData_x[i])

vSqdMean = np.mean(vSqd)

y = []
for i in range(0,3000):
    y.append(vSqdMean)
inTemp = np.mean(gasTempDataIn)
outTemp = np.mean(gasTempDataOut)

print("GasIn: " + str(inTemp))
print("GasOut " + str(outTemp))
print("T_COM: " + str(2./3. * vSqdMean))

# plt.plot(gasTempDataIn)
# plt.plot(gasTempDataOut)
# plt.plot(x,y)

plt.figure(1)
plt.plot(vSqd)
plt.plot((0,700),(vSqdMean,vSqdMean))
plt.show()




