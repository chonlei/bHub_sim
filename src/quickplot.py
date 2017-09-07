import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! - Otherwise figures try to load but -X not necessarily on
import matplotlib.pyplot as plt
import sys

fileN = sys.argv[1]
N = int(sys.argv[2]) + 1
S = np.loadtxt(fileN)
plt.plot(np.arange(N)/np.float(N)*100.,S)
plt.savefig("temp.png")
