# -*- coding: utf-8 -*-
"""
Created on Sun Oct 01 15:11:58 2017

@author: Jacky Lei
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
import os
import glob
import re

## figure setup
params = {
   'axes.labelsize': 14,
   'font.size': 14,
   'legend.fontsize': 10,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'text.usetex': False,
   'figure.figsize': [9, 6.6]
   }
plt.rcParams.update(params)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel(r"x [$\mu$m]")
ax.set_ylabel(r"y [$\mu$m]")
ax.set_zlabel(r"z [$\mu$m]")


CASE=101
iBatch = 67

###############################################################################
# setup data
try:
    fileName = sys.argv[1]
except Exception:
    fileName = None
try:
    isImagedCells = bool(int(sys.argv[2]))
except Exception:
    isImagedCells = True
try:
    nBatch = int(sys.argv[3])
except Exception:
    nBatch = 78 #len(glob.glob("case%s/*"%CASE)) - 3

if fileName == None:
    allfiles=glob.glob("../plot-traces/case%s/*"%CASE)
    fileName = [a for a in allfiles if "_p_1_" in a][0]
fileDir = os.path.dirname(fileName)
fileBase = os.path.basename(fileName)
fileIdx = re.findall(".*model_(\d+)_morphology_(\d+)_seed_(\d+)_mode_(\d+)_.*",fileName)
fileId = "model_%s_morphology_%s_seed_%s_mode_%s"%fileIdx[0]
mode = int(fileIdx[0][3])
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
varName = r"[Ca]$_i$ [$\mu$M]" if "Ca" in fileName else r"$V_m$ [mV]"
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

shapeX = shapeXBatch[iBatch]
X = np.memmap(fileNameBatch[iBatch], dtype='float64', mode='r', shape=shapeX)*1000.


###############################################################################
morphology = 5
pathToFile = "./"
coorDataFile = "H4-1-%d"%morphology
CoorData = np.loadtxt(pathToFile+coorDataFile+".txt")
# proceed CoorData to be acceptable format in genCoupleMatrix()
CoorData = CoorData[CoorData[:,0]==2][:,1:4]
centreCoorData = np.mean(CoorData,0)
CoorData = CoorData - centreCoorData

###############################################################################
xs,ys,zs = CoorData.T
values = X[:,0]

p = ax.scatter3D(-zs, ys, zs=xs, c=values, cmap='hot', s=30)

cb = fig.colorbar(p, ax=ax)
cb.set_label(r'[Ca]$_i$ [$\mu$M]')

plt.savefig("../figures/3dheatmap-%d.png"%CASE,bbox_inches='tight')
plt.savefig("../figures/3dheatmap-%d.pdf"%CASE,format='pdf',bbox_inches='tight')