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
varName = r"[Ca]$_i$ [mM]" if "Ca" in fileName else r"$V_m$ [mV]"
title = "short simulation homogeneous GJ @ mouse 40-3"
if isImagedCells:
    saveName = os.path.join(fileDir,'imaged_'+fileBase[:-4]+'.png')
else:
    saveName = os.path.join(fileDir,'whole_'+fileBase[:-4]+'.png')
XX = np.zeros((nBatch+1)*shapeX[1])

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
#fig = plt.figure()
#ax = fig.add_subplot(111)
aveX = np.zeros(X[0].shape)
maxCa = 0
for i in xrange(len(X)):
    if isImagedCells:
        if (i not in hubList) and (i in imagedCells):
            #ax.plot(t, X[i], 'k', alpha=0.3)
            aveX += X[i]
    else:
        if i not in hubList:
            #ax.plot(t, X[i], 'k', alpha=0.3)
            aveX += X[i]
#ax.plot(t, X[i], 'k', alpha=0.3, label='non-hub')
for i in np.array(hubList)[np.array(hubList)>0]:
    #ax.plot(t, X[i], 'r')
    aveX += X[i]
#ax.plot(t, X[i], 'r', label='hub')
aveX = aveX/float(shapeX[0]) if not isImagedCells else aveX/float(len(imagedCells))
#ax.plot(t, aveX, 'g', linewidth=3.0, label='average all')
maxCa = max(np.max(X),maxCa)
XX[0:len(aveX)] = aveX
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
                    #ax.plot(t, X[i], 'k', alpha=0.3)
                    aveX += X[i]
            else:
                if i not in hubList:
                    #ax.plot(t, X[i], 'k', alpha=0.3)
                    aveX += X[i]
        for i in np.array(hubList)[np.array(hubList)>0]:
            #ax.plot(t, X[i], 'r')
            aveX += X[i]
        #ax.plot(t, X[i], 'r')
        aveX = aveX/float(shapeX[0]) if not isImagedCells else aveX/float(len(imagedCells))
        #ax.plot(t, aveX, 'g', linewidth=3.0)
        maxCa = max(np.max(X),maxCa)
        #print len(XX[(iBatch+1)*len(aveX):(iBatch+2)*len(aveX)]), len(aveX), iBatch
        XX[(iBatch+1)*len(aveX):(iBatch+2)*len(aveX)] = aveX

# show where silencing is applied
if mode!=0:
    try:
        import matplotlib.patches as patches
        startName = "#silenceStart = "
        with open(os.path.join(fileDir,fileId+'.log'),"r") as fi:
            for ln in fi:
                if ln.startswith(startName):
                    temp = ln[len(startName):]
        silenceStart = float(temp)
        startName = "#silenceDur = "
        with open(os.path.join(fileDir,fileId+'.log'),"r") as fi:
            for ln in fi:
                if ln.startswith(startName):
                    temp = ln[len(startName):]
        silenceDur = float(temp)
        if maxCa>0.0005:
            rect_y = (0.00095,0.00001) if "Ca" in fileName else (-5, 1)
        else:
            rect_y = (0.000475,0.000005) if "Ca" in fileName else (-5, 1)
        #ax.add_patch(patches.Rectangle( (silenceStart,rect_y[0]), silenceDur, rect_y[1] , alpha=0.6))
    except Exception:
        pass

####################################
idxx=int(sys.argv[4])
stepsizehub=int(sys.argv[6])
#np.savetxt("tmp%d.txt"%idxx,XX)



sStartIdx = int(silenceStart/tstep)
sEndIdx = int((silenceStart+silenceDur)/tstep)
skipIdx = int(75*1000/tstep)
gStartIdx = int(50*1000/tstep)
#XX0 = XX[gStartIdx+skipIdx:sStartIdx-1]
XX1 = XX[sStartIdx+skipIdx:sEndIdx-1]
#minxx0 = np.min(XX0)
minxx1 = np.min(XX1)
#XX0 -= minxx0
XX1 -= minxx1
outpercent = np.mean(XX1)

inhibPer = idxx*stepsizehub/np.float(len(hubList))*100.
nameid = sys.argv[5]
with open("sim%s.txt"%nameid, 'a') as f:
    f.write("%f %f \n"%(inhibPer,outpercent))
print sStartIdx,sEndIdx,skipIdx,gStartIdx

# plotting setup
#ax.set_xlabel("t [ms]")
#ax.set_ylabel(varName)
#if "Ca" in fileName:
#    if maxCa>0.0005:
#        ax.set_ylim([0, 0.001])
#    else:
#        ax.set_ylim([0, 0.0005])
#else:
#    ax.set_ylim([-100.0, 0.0])
#ax.set_title(title)
#plt.savefig(saveName)

