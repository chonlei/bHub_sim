import numpy as np
import sys
import os
import modelSetup

"""
Args:
    fileName (str): path to log file
    isImagedCells (bool): 1 or 0; only analyse imaged cells
    isPrintOutHubs (bool): print out all the hubs defined in simulation
    isInteractivePlot (bool): use python interactive plotting interface
    cellAnalyse (int): index of cell to analyse
"""

# setup data
try:
    fileName = sys.argv[1]
except Exception:
    raise Exception("Please enter path to log file")

try:
    isImagedCells = bool(int(sys.argv[2]))
except Exception:
    isImagedCells = False

try:
    isPrintOutHubs = bool(int(sys.argv[3]))
except Exception:
    isPrintOutHubs = False

try:
    isInteractivePlot = bool(int(sys.argv[4]))
except Exception:
    isInteractivePlot = False

try:
    cellAnalyse = int(sys.argv[5])
except Exception:
    cellAnalyse = None  # We will just analyse the silenced cell


## same as getSpatialDistribution.py
if not isInteractivePlot:
    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! - Otherwise figures try to load but -X not necessarily on
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt


fileDir = os.path.dirname(fileName)
fileBase = os.path.basename(fileName)
if isImagedCells:
    saveName = os.path.join(fileDir,'imagedSD_'+fileBase[:-4]+'.png')
else:
    saveName = os.path.join(fileDir,'wholeSD_'+fileBase[:-4]+'.png')


# hubList (if have)
hubList = np.loadtxt(fileName,delimiter=',',dtype=int)

# get only imaged cells (if applicable)
try:
    imagedCells = []
    startName = "#imagedCells = "
    with open(fileName,"r") as fi:
        for ln in fi:
            if ln.startswith(startName):
                imagedCells = ln[len(startName):]
    imagedCells = [int(x.strip()) for x in imagedCells.split(',')]
    imagedHubs = []
    if isImagedCells:
        startName = "#imagedHubs = "
        with open(fileName,"r") as fi:
            for ln in fi:
                if ln.startswith(startName):
                    imagedHubs = ln[len(startName):]
        imagedHubs = [int(x.strip()) for x in imagedHubs.split(',')]
        #TODO: fix the index for those imaged cells only
        hubList = imagedHubs
except Exception:
    imagedCells = []
print hubList

startName = "#silencedCell = "
with open(fileName,"r") as fi:
    for ln in fi:
        if ln.startswith(startName):
            temp = ln[len(startName):]
silencedCell = [int(x.strip()) for x in temp.split(',')]

# get morphology info
startName = "#morphology = "
with open(fileName,"r") as fi: 
    for ln in fi: 
        if ln.startswith(startName):
            temp = ln[len(startName):]
morphology = int(temp)

startName = "#species = "
with open(fileName,"r") as fi:
    for ln in fi:
        if ln.startswith(startName):
            temp = ln[len(startName):]
species = int(temp)

startName = "#dthres = "
with open(fileName,"r") as fi:
    for ln in fi:
        if ln.startswith(startName):
            temp = ln[len(startName):]
dthres = float(temp)

startName = "#isletsize = "
with open(fileName,"r") as fi:
    for ln in fi:
        if ln.startswith(startName):
            temp = ln[len(startName):]
try:
    isletsize = float(temp)
except Exception:
    isletsize = None

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
if isPrintOutHubs:
    for i in hubList:
        print("cell %d ; links %d; coor %f,%f,%f"%(i,numLinks[i],CoorData[i,0],CoorData[i,1],CoorData[i,2]))

######
## Plot the spatial connectivity distribution
######
plt.figure(1)
countNumLinks = np.histogram(numLinks, bins=np.arange(numLinks.min(), numLinks.max()+1)-0.5)[0]
labelname = "islet" if not isImagedCells else "imaged" 
plt.plot(np.arange(numLinks.min(), numLinks.max()), countNumLinks, '-s',label=labelname)
plt.hist(numLinks[hubList],label='hubs')
plt.xlabel("#GJ to neighbouring cells")
plt.ylabel("#cells")
plt.legend()

#plt.savefig(saveName)

############################## until here #####################################
if cellAnalyse == None:
    cellAnalyse = silencedCell[0]
print("#### analysing cell %d"%cellAnalyse)
print("#GJ connection %d"%numLinks[cellAnalyse])
print("Connecting to cells: ")
CoupledMatrix = CoupledMatrix+CoupledMatrix.T
connectedCells = np.arange(ncells)[CoupledMatrix[:,cellAnalyse]>0]
print(connectedCells)
isConnectedCellsHub = [True if i in hubList else False for i in connectedCells] 
print(isConnectedCellsHub)
print("Its connecting cells connect to: ")
nextLayerConnectedCells = []
isNextLayerConnectedCellsHub = []
for i in connectedCells:
    temptemp = list(np.arange(ncells)[CoupledMatrix[:,i]>0])
    temp = list(temptemp)
    temp2 = []
    temp.remove(cellAnalyse)
    for j in temptemp:
        if j in connectedCells:
            temp.remove(j)
    for j in temp:
        temp2.append(True if j in hubList else False) 
    nextLayerConnectedCells.append(temp)
    isNextLayerConnectedCellsHub.append(temp2)
print(nextLayerConnectedCells)
print("Of which contains other hubs: ")
print(isNextLayerConnectedCellsHub)

# plot the islet
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
centreCoorData = np.mean(CoorData,0)
if isImagedCells:
    CoorDataMask = (np.sum((CoorData - centreCoorData)**2,1)<(isletsize**2))
else:
    CoorDataMask = (np.sum((CoorData - centreCoorData)**2,1)>0)
print "imagedCells: ", imagedCells
for i,(xs,ys,zs) in enumerate(CoorData[CoorDataMask,:]): #np.array(len(CoorData))[CoorDataMask]
    ax.scatter(xs, ys, zs, c='b', marker='o')
    if i in imagedCells:
        ax.scatter(xs, ys, zs, c='r', marker='o')
    if i == cellAnalyse:
        ax.scatter(xs, ys, zs, c='g', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# plot cells around cellAnalyse
secondLayer = [item for sublist in nextLayerConnectedCells for item in sublist]
fig = plt.figure(3)
ax = fig.add_subplot(111, projection='3d')
CoorDataMask = np.arange(ncells)[(np.sum((CoorData - CoorData[cellAnalyse,:])**2,1)<((3*dthres)**2))]
size = 60
for i,(xs,ys,zs) in enumerate(CoorData[:,:]):
    if i in CoorDataMask:
        ax.scatter(xs, ys, zs, c='b', marker='o')
        if i == cellAnalyse:
            ax.scatter(xs, ys, zs, c='g', marker='o',s=size)
        if i in connectedCells:
            ax.scatter(xs, ys, zs, c='c', marker='o',s=size)
        if i in secondLayer:
            ax.scatter(xs, ys, zs, c='y', marker='o',s=size)
    if i in hubList and i!=cellAnalyse:
        ax.scatter(xs, ys, zs, c='r', marker='o',s=size)
plt.savefig("temp.png")
if isInteractivePlot:
    plt.show()
