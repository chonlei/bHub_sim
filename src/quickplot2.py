import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! - Otherwise figures try to load but -X not necessarily on
import matplotlib.pyplot as plt
import sys
import os
import re

fileName = sys.argv[3]
fileDir = os.path.dirname(fileName)
fileBase = os.path.basename(fileName)
fileIdx = re.findall(".*model_(\d+)_morphology_(\d+)_seed_(\d+)_mode_(\d+)_.*",fileName)
fileId = "model_%s_morphology_%s_seed_%s_mode_%s"%fileIdx[0]

fileN = sys.argv[1]
fileID = sys.argv[4]
N = int(sys.argv[2]) + 1
try:
        hubList = np.loadtxt(os.path.join(fileDir,fileId+'.log'),delimiter=',',dtype=int)
except Exception:
        hubList = []
S = np.loadtxt(fileN)
#plt.plot(np.arange(N)*6./np.float(len(hubList))*100.,S/S[0])
plt.plot(S[:,0],S[:,1]/S[0,1])
plt.savefig("sim%s.png"%fileID)
