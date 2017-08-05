import re
import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! - Otherwise figures try to load but -X not necessarily on
import matplotlib.pyplot as plt
import sys
import os
import glob

# setup data
try:
    fileName = sys.argv[1]
except Exception:
    fileName = "Ca_WT_5phubs_9s_sGJ_100by501.dat"
try:
    isImagedCells = bool(int(sys.argv[2]))
except Exception:
    isImagedCells = False
try:
    nBatch = int(sys.argv[3])
except Exception:
    nBatch = 0
fileDir = os.path.dirname(fileName)
fileBase = os.path.basename(fileName)
fileIdx = re.findall(".*model_(\d+)_morphology_(\d+)_seed_(\d+)_mode_(\d+)_.*",fileName)
fileId = "model_%s_morphology_%s_seed_%s_mode_%s"%fileIdx[0]
shape = re.findall('.*_(\w+)x(\w+)\.dat',fileName)[0]
shapeX = (int(shape[0]),int(shape[1]))
if nBatch>0:
    filePrefix,startBatch = re.findall("(.*)_p_(\d+)_.*",fileName)[0]
    fileNameBatch = []
    shapeXBatch = []
    for i in range(1,nBatch+1):
        getBatch = int(startBatch)+i
        tempfilename = glob.glob(filePrefix+"_p_%d_*"%getBatch)[0]
        tempshape = re.findall('.*_(\w+)x(\w+)\.dat',tempfilename)[0]
        shapeXBatch.append((int(tempshape[0]),int(tempshape[1])))
        fileNameBatch.append(tempfilename)
varName = r"[Ca]$_i$ [mM]" if "Ca" in fileName else r"$V_m$ [mV]"
title = "short simulation homogeneous GJ @ mouse 40-3"
if isImagedCells:
    saveName = os.path.join(fileDir,'imaged_'+fileBase[:-4]+'.png')
else:
    saveName = os.path.join(fileDir,'whole_'+fileBase[:-4]+'.png')

# load data
startName = "#dt = "
with open(os.path.join(fileDir,fileId+'.log'),"r") as fi:
    for ln in fi:
        if ln.startswith(startName):
            temp = ln[len(startName):]
dt = float(temp)
startName = "#downSampling = "
with open(os.path.join(fileDir,fileId+'.log'),"r") as fi:
    for ln in fi:
        if ln.startswith(startName):
            temp = ln[len(startName):]
downSampling = float(temp)
tstep = dt*downSampling
X = np.memmap(fileName, dtype='float64', mode='r', shape=shapeX)
t = np.arange(0,shapeX[1]*tstep,tstep)  # account dt and downSampling

# hubList (if have)
try:
    hubList = np.loadtxt(os.path.join(fileDir,fileId+'.log'),delimiter=',',dtype=int)
except Exception:
    hubList = []

# get only imaged cells (if applicable)
if isImagedCells:
    imagedCells = []
    startName = "#imagedCells = "
    with open(os.path.join(fileDir,fileId+'.log'),"r") as fi:
        for ln in fi:
            if ln.startswith(startName):
                imagedCells = ln[len(startName):]
    imagedCells = [int(x.strip()) for x in imagedCells.split(',')]
    imagedHubs = []
    startName = "#imagedHubs = "
    with open(os.path.join(fileDir,fileId+'.log'),"r") as fi:
        for ln in fi:
            if ln.startswith(startName):
                imagedHubs = ln[len(startName):]
    imagedHubs = [int(x.strip()) for x in imagedHubs.split(',')]
    hubList = imagedHubs

# main plot
aveX = np.zeros(X[0].shape)
for i in xrange(len(X)):
    if isImagedCells:
        if (i not in hubList) and (i in imagedCells):
            plt.plot(t, X[i], 'k', alpha=0.3)
            aveX += X[i]
    else:
        if i not in hubList:
            plt.plot(t, X[i], 'k', alpha=0.3)
            aveX += X[i]
plt.plot(t, X[i], 'k', alpha=0.3, label='non-hub')
for i in np.array(hubList)[np.array(hubList)>0]:
    plt.plot(t, X[i], 'r')
    aveX += X[i]
plt.plot(t, X[i], 'r', label='hub')
aveX = aveX/float(shapeX[0]) if not isImagedCells else aveX/float(len(imagedCells))
plt.plot(t, aveX, 'g', linewidth=3.0, label='average all')
# plot the rest if splitted into batches
if nBatch>0:
    for iBatch in range(nBatch):
        shapeX = shapeXBatch[iBatch]
        X = np.memmap(fileNameBatch[iBatch], dtype='float64', mode='r', shape=shapeX)
        t = np.arange(t[-1],t[-1]+shapeX[1]*tstep,tstep)  # account dt
        aveX = np.zeros(X[0].shape)
        for i in xrange(len(X)):
            if isImagedCells:
                if (i not in hubList) and (i in imagedCells):
                    plt.plot(t, X[i], 'k', alpha=0.3)
                    aveX += X[i]
            else:
                if i not in hubList:
                    plt.plot(t, X[i], 'k', alpha=0.3)
                    aveX += X[i]
        for i in np.array(hubList)[np.array(hubList)>0]:
            plt.plot(t, X[i], 'r')
            aveX += X[i]
        plt.plot(t, X[i], 'r')
        aveX = aveX/float(shapeX[0]) if not isImagedCells else aveX/float(len(imagedCells))
        plt.plot(t, aveX, 'g', linewidth=3.0)


# plotting setup
plt.xlabel("t [ms]")
plt.ylabel(varName)
if "Ca" in fileName:
    plt.ylim([0, 0.0005])
#plt.title(title)
plt.savefig(saveName)

