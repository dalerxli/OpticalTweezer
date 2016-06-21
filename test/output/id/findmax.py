import numpy as np
import scipy as sp
import re
import matplotlib.pyplot as plot

def findmaxvalue(number):
    if(number < 10):
        filename = "100"+str(number)+".dat"
    elif(number >= 10 and number < 100):
        filename = "10"+str(number)+".dat"
    elif(number >= 100):
        filename = "1"+str(number)+".dat"

    data_array = []
    with open(filename,"r") as infile:
        line = infile.read()
        data = line.split("\n")
        for i in data:
            if i != "test" and i != "1":
                data_array.append(i)

    for i in range(0,len(data_array)):
        data_array[i] = data_array[i].split("\t")


    for i in range(0,len(data_array)-1):
        data_array[i].remove("GasIn")

    data_array.remove(data_array[len(data_array)-1])


    for i in range(0,len(data_array)):
        for j in range(0,len(data_array[i])):
            data_array[i][j] = float(data_array[i][j])
    x_data  = []
    y_data  = []
    z_data  = []
    vx_data = []
    vy_data = []
    vz_data = []
    v_data  = []
    ax_data = []
    ay_data = []
    az_data = []
    a_data  = []

    for i in range(0,len(data_array)):
        x_data.append(data_array[i][0])
        y_data.append(data_array[i][1])
        z_data.append(data_array[i][2])
        vx_data.append(data_array[i][3])
        vy_data.append(data_array[i][4])
        vz_data.append(data_array[i][5])
        v_data.append(data_array[i][6])
        ax_data.append(data_array[i][7])
        ay_data.append(data_array[i][8])
        az_data.append(data_array[i][9])
        a_data.append(data_array[i][10])

    max_a= np.amax(a_data)
    max_a_index = np.argmax(a_data)
    return([max_a_index,max_a])


vals = []
for i in range(0,500):
    vals.append(findmaxvalue(i))
for i in range(0,len(vals)):
    print(str(vals[i][0]) + "   " + str(vals[i][1]))
# print(vals)
# plot.scatter(vals)
# plot.show()
