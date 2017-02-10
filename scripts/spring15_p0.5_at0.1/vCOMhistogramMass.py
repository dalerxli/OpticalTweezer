import math
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from subprocess import call
from scipy.stats import norm

# proc = call("ls *.dat",shell=True)
# datetime = "170123_2033_"
datetime = sys.argv[1]+"_"
gasTempDataIn = np.genfromtxt(datetime+"gasTempData.dat",usecols=0,skip_header=100)
gasTempDataOut = np.genfromtxt(datetime+"gasTempData.dat",usecols=1,skip_header=100)
vCOMData_x = np.genfromtxt(datetime+"vCOMData.dat",usecols=0,skip_header=100)
vCOMData_y = np.genfromtxt(datetime+"vCOMData.dat",usecols=1,skip_header=100)
vCOMData_z = np.genfromtxt(datetime+"vCOMData.dat",usecols=2,skip_header=100)
internalTempData = np.genfromtxt(datetime+"temperature_internal.dat",skip_header=200)
N = 32

internalTemp = np.mean(internalTempData)

vSqd = []
for i in range(0,len(vCOMData_x)):
    vSqd.append(32*(vCOMData_x[i]*vCOMData_x[i]+vCOMData_x[i]*vCOMData_x[i]+vCOMData_x[i]*vCOMData_x[i])*0.5)

vSqdMean = np.mean(vSqd)

# histogram_x,bins_x = np.histogram(vCOMData_x,bins=100,normed=False)
# histogram_y,bins_y = np.histogram(vCOMData_y,bins=100,normed=False)
# histogram_z,bins_z = np.histogram(vCOMData_z,bins=100,normed=False)
histogram_x,bins_x = np.histogram(vCOMData_x,bins=100,normed=True)
histogram_y,bins_y = np.histogram(vCOMData_y,bins=100,normed=True)
histogram_z,bins_z = np.histogram(vCOMData_z,bins=100,normed=True)

inTemp = np.mean(gasTempDataIn)
outTemp = np.mean(gasTempDataOut)

statistics = open(datetime+"statistics_mass.dat","w")
statistics.write("GasIn: " + str(inTemp)+"\n")
statistics.write("GasOut: " + str(outTemp)+"\n")
statistics.write("T_COM: " + str(2./3. * vSqdMean)+" +- " + str(np.std(vSqd)) + "\n")
statistics.write("T_INT: " + str(internalTemp)+"\n")

statistics.write("Mu_x " + str(np.mean(vCOMData_x))+"\n")
statistics.write("Sigma_x: " + str(np.std(vCOMData_x))+"\n")
statistics.write("Mu_y " + str(np.mean(vCOMData_y))+"\n")
statistics.write("Sigma_y: " + str(np.std(vCOMData_y))+"\n")
statistics.write("Mu_z " + str(np.mean(vCOMData_z))+"\n")
statistics.write("Sigma_z: " + str(np.std(vCOMData_z))+"\n")

histogram_x_file = open(datetime+"histogram_mass_vx.dat","w")
histogram_y_file = open(datetime+"histogram_mass_vy.dat","w")
histogram_z_file = open(datetime+"histogram_mass_vz.dat","w")
for i in range(0,len(histogram_x)):
    histogram_x_file.write(str(bins_x[i]) + "\t" + str(histogram_x[i]) + "\n")
    histogram_y_file.write(str(bins_y[i]) + "\t" + str(histogram_y[i]) + "\n")
    histogram_z_file.write(str(bins_z[i]) + "\t" + str(histogram_z[i]) + "\n")


# plt.figure(1)
# plt.hist(vCOMData_x,bins=100)
# plt.figure(2)
# plt.hist(vCOMData_y,bins=100)
# plt.figure(3)
# plt.hist(vCOMData_z,bins=100)
# plt.show()
# plt.figure(1)
# plt.plot(vSqd)
# plt.plot((0,700),(vSqdMean,vSqdMean))
# plt.figure(2)
# plt.hist(vCOMData_x,bins=100,normed=True)
# plt.plot(x,gasInPDF)
# plt.show()
