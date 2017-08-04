import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! - Otherwise figures try to load but -X not necessarily on
import matplotlib.pyplot as plt
import sys
import os
import modelSetup


# setup data
try:
    fileName = sys.argv[1]
except Exception:
    raise Exception("Please enter path to log file")

try:
    isImagedCells = bool(sys.argv[2])
except Exception:
    isImagedCells = False

fileDir = os.path.dirname(fileName)
fileBase = os.path.basename(fileName)
if isImagedCells:
    saveName = os.path.join(fileDir,'imagedSD_'+fileBase[:-4]+'.png')
else:
    saveName = os.path.join(fileDir,'wholeSD_'+fileBase[:-4]+'.png')


# hubList (if have)
hubList = np.loadtxt(fileName,delimiter=',',dtype=int)

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

startName = "#silencedCell = "
with open(os.path.join(fileDir,fileId+'.log'),"r") as fi:
    for ln in fi:
        if ln.startswith(startName):
            temp = ln[len(startName):]
silencedCell = [int(x.strip()) for x in temp.split(',')]

# get morphology info
startName = "#morphology = "
with open(os.path.join(fileDir,fileId+'.log'),"r") as fi: 
    for ln in fi: 
        if ln.startswith(startName):
            temp = ln[len(startName):]
morphology = int(temp)

startName = "#species = "
with open(os.path.join(fileDir,fileId+'.log'),"r") as fi:
    for ln in fi:
        if ln.startswith(startName):
            temp = ln[len(startName):]
species = int(temp)

if species==0:
    #pathToCoupledMatrix = '../morphologies/mouse/CouplingMatrix-mouse40-3-175.dat'
    #pathToCoupledMatrix = '../morphologies/mouse/CouplingMatrixMouse403.dat'
    pathToMorphology = "../morphologies/mouse/Mouse 40-%d.txt"%morphology
elif species==1:
    pathToMorphology = ""
elif species==2:
    pathToMorphology = "../morphologies/cubic_lattice/cubic%d.txt"%morphology


######
## Import system setup files (.hoc files and system matrix)
######
#CoupledMatrix = np.loadtxt(pathToCoupledMatrix,delimiter=' ')#[-100:,-100:] #only taking last 100 cells
CoorData = np.loadtxt(pathToMorphology)
# process CoorData to be acceptable format in modelSetup.genCoupleMatrix()
CoorData = CoorData[CoorData[:,0]==11][:,1:4]
CoupledMatrix = modelSetup.genCoupleMatrix(CoorData,dthres,isletsize,True)
ncells = CoupledMatrix.shape[0]
if isImagedCells:
    CoupledMatrix = CoupledMatrix[imagedCells,:][:,imagedCells]
numLinks = np.sum(CoupledMatrix+CoupledMatrix.T,1)


######
## Plot the spatial connectivity distribution
######
plt.figure(1)
countNumLinks = np.histogram(numLinks, bins=np.arange(numLinks.min(), numLinks.max()+1)-0.5)[0]
labelname = "islet" if not isImagedCells else "imaged" 
plt.plot(np.arange(numLinks.min(), numLinks.max()), countNumLinks, '-s',label=labelname)
plt.hist(numLinks[hubList])
plt.xlabel("#GJ to neighbouring cells")
plt.ylabel("#cells")

plt.savefig(saveName)

