import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys

# import vCOMdata.dat as array
# folder="../output/runs/170123_2033/"
folder="../output/runs/"+str(sys.argv[1])
# when = str(sys.argv[1])
for i in range(0,len(sys.argv)):
    print(str(i) + ": "+ str(sys.argv[i]))
# data = np.genfromtxt(folder+"/vCOMData.dat",usecols=(0,1,2), skip_header=100)
# vx = np.genfromtxt(folder+" "+str(sys.argv[1])+"/vCOMData.dat",usecols=0, skip_header=100)
vx = np.genfromtxt(folder+"/vCOMData.dat",usecols=0, skip_header=100)
vy = np.genfromtxt(folder+"/vCOMData.dat",usecols=1, skip_header=100)
vz = np.genfromtxt(folder+"/vCOMData.dat",usecols=2, skip_header=100)
# print(vx)

# vx = []
# for i in range(0,len(data)):
    # vx.append(data[i][0])
# hist_vx = np.histogram(vx,bins=100)
# print(len(hist_vx[1]))
# print("")
# print(len(hist_vx[0]))
# hist_vy = np.histogram(vy,bins=100)
# hist_vz = np.histogram(vz,bins=100)

# printable = []
# for i in range(0,len(hist_vx[0])):
    # printable.append((hist_vx[1][i],hist_vx[0][i]))

# print(hist_vx)

plt.hist(vx,bins=100)
# plt.hist(vy,bins=100)
# plt.hist(vz,bins=100)
# plt.plot(hist_vx[1],hist_vx[0])
plt.show()
