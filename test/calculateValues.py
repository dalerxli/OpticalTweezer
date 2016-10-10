import numpy as np
import os

os.chdir("/home/motz/codes/OpticalTweezer/output/runs/160919_1756/gas")
files = os.listdir("./")
files_sorted = sorted(files)
for i in range(0,10):
    print(files[i])

values = []
for i in range(0,10000):
    open_file = open(files[i],"r")
    vals = np.loadtxt(open_file,usecols=(3,))
    if(len(vals) != 0):
        print(vals[0])
        values.append(vals[0])
    
    # values.append(vals[0])
    # values.append(np.genfromtxt(open_file,usecols=(3,),max_rows=2))

print(values)
print(np.mean(values))
print(np.var(values))
print(np.std(values))
print(len(values))
# open_file = open(files[1],"r")
# values = open_file.readline()
# values = values.split("\t")
# vals = np.loadtxt(open_file,usecols=(3,))
# print(vals[0])


