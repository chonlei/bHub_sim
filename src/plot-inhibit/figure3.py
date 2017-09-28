import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

try:
    CASE = int(sys.argv[1])
except Exception:
    raise Exception("Enter case idx to plot as arg")
SUBCASE = ['hub','follower']
colour = ['r','k']

for (idx,subcase) in enumerate(SUBCASE):
    loc = "case%d/%s/"%(CASE,subcase)  # location of files
    FILES = glob.glob("%s*.txt"%loc)  # list of files to plot
    for sim in FILES:
        trace = np.loadtxt(sim)
        trace[:,1] = trace[:,1]/trace[0,1]*100.0  # normalise
        plt.plot(trace[:,0],trace[:,1],c=colour[idx])

plt.xlabel("%hub inhibited")
plt.ylabel("%activity")
plt.show()
