import math
import numpy as np
import re

outputfile = open("gastemp_plot.dat","w")
folderlistdata = open("folderlist.txt","r")
folderlist = folderlistdata.read()
folders = folderlist.split("\n")
# paramtest = open(folders[1]+"_parameters.txt","r")
# paramdata = ""
# for line in paramtest:
    # if re.search("dQ",line):
        # paramdata = line

# dQdata = paramdata.split()
# print dQdata[2]

for i in range(0,len(folders)-1):
    parameterdata = open(folders[i]+"_parameters.txt")
    for line in parameterdata:
        if re.search("dQ",line):
            paramdata = line
    dQdata = paramdata.split()
    gasInData = np.genfromtxt(folders[i]+"_gasTempData.dat",usecols=0,skip_header=100)
    gasOutData = np.genfromtxt(folders[i]+"_gasTempData.dat",usecols=1,skip_header=100)
    inTemp = np.mean(gasInData)
    outTemp = np.mean(gasOutData)
    tcm = (math.pow(inTemp,3./2.)+math.pi/8.*math.pow(outTemp,3./2.))/(math.sqrt(inTemp)+math.pi/8.*math.sqrt(outTemp))
    outputfile.write(str(dQdata[2] + "\t" + str(tcm) + "\n"))

