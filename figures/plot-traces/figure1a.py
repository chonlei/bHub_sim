import re
import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! - Otherwise figures try to load but -X not necessarily on
import matplotlib.pyplot as plt
import sys
import os
import glob

CASE=6


## figure setup
params = {
   'axes.labelsize': 12,
   'font.size': 12,
   'legend.fontsize': 10,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'text.usetex': False,
   'figure.figsize': [6.5, 4.5]
   }
plt.rcParams.update(params)

#fig = plt.figure(figsize=(13.5, 4.5))
#ax = fig.add_subplot(1, 1, 1)
f, (ax, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
f.subplots_adjust(hspace=0)

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
    nBatch = 40 #len(glob.glob("case%s/*"%CASE)) - 3

if fileName == None:
    allfiles=glob.glob("case%s/*"%CASE)
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
X = np.memmap(fileName, dtype='float64', mode='r', shape=shapeX)*1000.
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
maxCa = 0
for i in xrange(len(X)):
    if isImagedCells:
        if (i not in hubList) and (i in imagedCells):
            #ax.plot(t/1000., X[i], 'k', alpha=0.3)
            aveX += X[i]
    else:
        if i not in hubList:
            ax.plot(t/1000., X[i], 'k', alpha=0.3)
            aveX += X[i]
#ax.plot(t/1000., X[i], 'k', alpha=0.3, label='non-hub')
for i in np.array(hubList)[np.array(hubList)>0]:
    if i == 1:
        ax2.plot(t/1000., X[i], 'r')
    aveX += X[i]
#ax.plot(t/1000., X[i], 'r', label='hub')
aveX = aveX/float(shapeX[0]) if not isImagedCells else aveX/float(len(imagedCells))
ax3.plot(t/1000., aveX, 'g', linewidth=3.0, label='average all')
maxCa = max(np.max(X),maxCa)
# plot the rest if splitted into batches
Xall = np.zeros((len(imagedCells),shapeXBatch[0][1]*nBatch))
tall = np.zeros(shapeXBatch[0][1]*nBatch)
if nBatch>0:
    for iBatch in range(nBatch):
        counter = 0
        shapeX = shapeXBatch[iBatch]
        X = np.memmap(fileNameBatch[iBatch], dtype='float64', mode='r', shape=shapeX)*1000.
        t = np.arange(t[-1],t[-1]+shapeX[1]*tstep,tstep)  # account dt
        tall[iBatch*shapeX[1]:(iBatch+1)*shapeX[1]] = t
        aveX = np.zeros(X[0].shape)
        for i in xrange(len(X)):
            if isImagedCells:
                if (i not in hubList) and (i in imagedCells):
                    #ax.plot(t/1000., X[i], 'k', alpha=0.3)
                    aveX += X[i]
                    Xall[counter,iBatch*shapeX[1]:(iBatch+1)*shapeX[1]] = X[i]
                    counter+=1
            else:
                if i not in hubList:
                    ax.plot(t/1000., X[i], 'k', alpha=0.3)
                    aveX += X[i]
        for i in np.array(hubList)[np.array(hubList)>0]:
            if i == 1:
                ax2.plot(t/1000., X[i], 'r')
            Xall[counter,iBatch*shapeX[1]:(iBatch+1)*shapeX[1]] = X[i]
            counter+=1
            aveX += X[i]
        #ax.plot(t/1000., X[i], 'r')
        aveX = aveX/float(shapeX[0]) if not isImagedCells else aveX/float(len(imagedCells))
        ax3.plot(t/1000., aveX, 'g', linewidth=3.0)
        maxCa = max(np.max(X),maxCa)
ax3.plot(tall/1000., np.mean(Xall,0)+np.std(Xall,0), 'g', linewidth=1.0)
ax3.plot(tall/1000., np.mean(Xall,0)-np.std(Xall,0), 'g', linewidth=1.0)
maxCa=0






###############################################################################
if True:
    allfiles=glob.glob("case%s/*"%CASE)
    fileName = [a for a in allfiles if "_p_0_" in a][1]
fileDir = os.path.dirname(fileName)
fileBase = os.path.basename(fileName)
fileIdx = re.findall(".*model_(\d+)_morphology_(\d+)_seed_(\d+)_mode_(\d+)_.*",fileName)
fileId = "model_%s_morphology_%s_seed_%s_mode_%s"%fileIdx[0]
mode = int(fileIdx[0][3])
shape = re.findall('.*_(\w+)x(\w+)\.dat',fileName)[0]
shapeX = (int(shape[0]),int(shape[1]))
nBatch += 1
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
#varName = r"[Ca]$_i$ [$\mu$M]" if "Ca" in fileName else r"$V_m$ [mV]"

# main plot
t = [0]
if nBatch>0:
    for iBatch in range(nBatch):
        shapeX = shapeXBatch[iBatch]
        X = np.memmap(fileNameBatch[iBatch], dtype='float64', mode='r', shape=shapeX)
        t = np.arange(t[-1],t[-1]+shapeX[1]*tstep,tstep)  # account dt
        ax.plot(t/1000., X[1], 'tab:blue')



###############################################################################
# show where silencing is applied
if mode!=0:
    try:
        import matplotlib.patches as patches
        startName = "#silenceStart = "
        with open(os.path.join(fileDir,fileId+'.log'),"r") as fi:
            for ln in fi:
                if ln.startswith(startName):
                    temp = ln[len(startName):]
        #silenceStart = float(temp)
        silenceStart = 45
        startName = "#silenceDur = "
        with open(os.path.join(fileDir,fileId+'.log'),"r") as fi:
            for ln in fi:
                if ln.startswith(startName):
                    temp = ln[len(startName):]
        silenceDur = 400 #float(temp)
        if maxCa>0.0005:
            rect_y = (0.95,0.01) if "Ca" in fileName else (-15, 10)
        else:
            rect_y = (0.525,0.025) if "Ca" in fileName else (-5.5, 3.5)
        ax.add_patch(patches.Rectangle( (silenceStart,rect_y[0]), silenceDur, rect_y[1] , alpha=0.6))
    except Exception:
        pass



# plotting setup
ax3.set_xlabel("t [s]")
ax.set_ylabel(r"$V_m$ [mV]")
ax2.set_ylabel(r"[Ca]$_i$ [$\mu$M]")
ax3.set_ylabel(r"E([Ca]$_i$) [$\mu$M]")

ax.set_ylim([-75.0, 2.5])
ax2.set_ylim([0.05, 0.5])
ax3.set_ylim([0.05, 0.5])
ax.set_xlim([0.0, 200])
'''
# add legend
legend = plt.legend()  #plt.legend(["case1", "case2", "case4"])#, loc=4);
frame = legend.get_frame()
frame.set_facecolor('0.9')
frame.set_edgecolor('0.9')
'''
plt.savefig("../figures/trace-%d_0.png"%CASE,bbox_inches='tight')
plt.savefig("../figures/trace-%d_0.pdf"%CASE,format='pdf',bbox_inches='tight')

plt.show()

