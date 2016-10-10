import numpy as np

data = np.loadtxt("first_lines.txt",usecols=(3,))

print(str(np.mean(data)) + "+- "+ str(np.std(data)))
print(np.max(data))
print(np.min(data))
