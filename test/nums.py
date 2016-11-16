import numpy as np

data = np.loadtxt("test.dat",dtype=int)
m = []
for j in range(0,30):
    m.append([])
print(m)
for i in data:
    m[i].append(int(i))

for j in range(0,30):
    print(len(m[j]))
