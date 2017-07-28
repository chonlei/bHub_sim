import re
import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! - Otherwise figures try to load but -X not necessarily on
import matplotlib.pyplot as plt
import sys
import os

# setup data
try:
    fileName = sys.argv[1]
except Exception:
    fileName = "Ca_WT_5phubs_9s_sGJ_100by501.dat"
fileDir = os.path.dirname(fileName)
fileIdx = re.findall(".*model_(\d+)_morphology_(\d+)_seed_(\d+)_mode_(\d+)_.*",fileName)
fileId = "model_%s_morphology_%s_seed_%s_mode_%s"%fileIdx[0]
shape = re.findall('.*_(\w+)x(\w+)\.dat',fileName)[0]
shapeX = (int(shape[0]),int(shape[1]))
varName = r"[Ca]$_i$ [mM]"
title = "short simulation homogeneous GJ @ mouse 40-3"
saveName = fileName[:-4]+'.png'

# load data
X = np.memmap(fileName, dtype='float64', mode='r', shape=shapeX)
time = np.arange(shapeX[1])*10.0  # account dt

# hubList (if have)
hubList = np.loadtxt(os.path.join(fileDir,fileId+'.log'),delimiter=',',dtype=int)

# main plot
aveX = np.zeros(X[0].shape)
for i in xrange(len(X)):
    if i not in hubList:
        plt.plot(time, X[i], 'k', alpha=0.3)
        aveX += X[i]
plt.plot(time, X[i], 'k', alpha=0.3, label='non-hub')
for i in np.array(hubList)[np.array(hubList)>0]:
    plt.plot(time, X[i], 'r')
    aveX += X[i]
plt.plot(time, X[i], 'r', label='hub')
aveX = aveX/float(shapeX[0])
plt.plot(time, aveX, 'g', linewidth=3.0, label='average all')

# plotting setup
plt.xlabel("t [ms]")
plt.ylabel(varName)
plt.title(title)
plt.savefig(saveName)

