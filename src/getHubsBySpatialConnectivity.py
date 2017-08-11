# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 23:24:03 2017

@author: Chon Lei
"""

import numpy as np
import matplotlib.pyplot as plt
import modelSetup
import sys
import warnings
import random

try:
    speices = int(sys.argv[1])
    morphology = int(sys.argv[2])
except Exception:
    raise Exception("Please enter speices ID and morphology ID")
try:
    meanGJ = int(sys.argv[3])
except Exception:
    meanGJ = 7
try:
    seed = int(sys.argv[4])
except Exception:
    seed = 1
try:
    dthres = float(sys.argv[5])
except Exception:
    dthres = 17.5
try:
    dcentre = np.float(sys.argv[6])
    isImitateExp = int(sys.argv[7])
except Exception:
    dcentre = 100000
    isImitateExp = 1

# Define data file and parameters for testing
if speices==0:
    pathToFile = "../morphologies/mouse/"
    coorDataFile = "Mouse 40-%d"%morphology
elif speices==1:
    #get human islet data
    pass
elif speices==2:
    pathToFile = "../morphologies/cubic_lattice/"
    coorDataFile = "cubic%d"%morphology

pHubs = 0.01
isInteractive = not False
rangeGJ = 1  # allows +/- rangeGJ of meanGJ
N_GJ = range(meanGJ-rangeGJ, meanGJ+rangeGJ)  # n GJ links that hubs allow to have
CoorData = np.loadtxt(pathToFile+coorDataFile+".txt")

random.seed(seed)
np.random.seed(seed)


def genCoupleMatrix(coorData, cutoff=20):
    # Generate a binary CoupleMatrix using islet coordinates data
    # 
    # if given 3 columns coordinate data of one cell type (e.g. b-cell)
    if isinstance(coorData, (str, unicode)):
        coorData = np.loadtxt(coorData)
    # else assume it is a numpy array type with equivalent format
    # expecting format of N-by-3 array (for N cells)
    N = np.shape(coorData)[0]        # change here if different format
    # create empty CoupleMatrix
    CoupleMatrix = np.zeros([N,N])
    for i in xrange(N):
        for j in xrange(N):
            diff = coorData[i,:] - coorData[j,:]
            CoupleMatrix[i,j] = np.inner(diff,diff)
    cutoff = cutoff**2            # compare in distance-square
    # convert it to binary matrix
    CoupleMatrix[CoupleMatrix<cutoff] = 1
    CoupleMatrix[CoupleMatrix>cutoff] = 0
    return CoupleMatrix


# proceed CoorData to be acceptable format in genCoupleMatrix()
CoorData = CoorData[CoorData[:,0]==11][:,1:4]
imagedCells = modelSetup.getImagedCellIdx(CoorData,topDir=2,imageDepth=10,Ncells=100,method=0)
cm = genCoupleMatrix(CoorData,dthres)
numLinks = np.sum(cm,1)
if isInteractive:
    plt.figure(1)
    countNumLinks = np.histogram(numLinks, bins=np.arange(numLinks.min(), numLinks.max()+1)-0.5)[0]
    plt.plot(np.arange(numLinks.min(), numLinks.max()), countNumLinks, '-s',label="$R_{th} = %g$"%dthres)
    
    plt.legend()
    plt.title("GJ spatial network distribution - %s"%coorDataFile)
    plt.xlabel("number of GJ to neighbouring cells")
    plt.ylabel("number of cells")
    #plt.savefig("SN_dist-%s.png"%coorDataFile)
    
    # turn cm to usable coupling matrix in NEURON script
    #cm = np.tril(cm,-1)
    #np.savetxt("CouplingMatrix-%s-%d.dat"%(coorDataFile,int(i*10),cm))
if isInteractive:
    plt.figure(3)
    cm1 = cm[imagedCells,:]
    numLinks1 = np.sum(cm1,1)
    countNumLinks1 = np.histogram(numLinks1, bins=np.arange(numLinks1.min(), numLinks1.max()+1)-0.5)[0]
    plt.plot(np.arange(numLinks1.min(), numLinks1.max()), countNumLinks1, '-s',label="all links")
    cm2 = cm[imagedCells,:][:,imagedCells]
    numLinks2 = np.sum(cm2,1)
    countNumLinks2 = np.histogram(numLinks2, bins=np.arange(numLinks2.min(), numLinks2.max()+1)-0.5)[0]
    plt.plot(np.arange(numLinks2.min(), numLinks2.max()), countNumLinks2, '-s',label="links amount imaged cells")
    plt.legend()
    plt.title("imaged cells GJ spatial network distribution - %s"%coorDataFile)
    plt.xlabel("number of GJ to neighbouring cells")
    plt.ylabel("number of cells")

# plot the islet
if isInteractive:
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure(2)
    ax = fig.add_subplot(111, projection='3d')
    centreCoorData = np.mean(CoorData,0)
    CoorDataMask = (np.sum((CoorData - centreCoorData)**2,1)<(dcentre**2))
    print "imagedCells: ", imagedCells
    for i,(xs,ys,zs) in enumerate(CoorData[CoorDataMask,:]): #np.array(len(CoorData))[CoorDataMask]
        ax.scatter(xs, ys, zs, c='b', marker='o')
        if i in imagedCells:
            ax.scatter(xs, ys, zs, c='r', marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

temp = []
numHubs = int(pHubs*len(numLinks))
for nGJ in N_GJ:
    temptemp = np.where(numLinks==nGJ)[0]
    temp += list(temptemp)
random.shuffle(temp)
hubsList = temp[0:numHubs]
imagedHubs = list(set(hubsList).intersection(imagedCells))
if len(imagedHubs) < int(pHubs*len(imagedCells)):
    nMorehubs = int(pHubs*len(imagedCells)) - len(imagedHubs)
    tempCellsToPick = [x for x in imagedCells if x not in imagedHubs]
    random.shuffle(tempCellsToPick)
    hubsList += tempCellsToPick[0:nMorehubs]
    imagedHubs += tempCellsToPick[0:nMorehubs]
print("hubsList: ",hubsList)
print("imagedHubs: ",imagedHubs)
for i in hubsList:
    print("cell %d; links %d"%(i,numLinks[i]))

if isInteractive:
    fig = plt.figure(4)
    ax = fig.add_subplot(111, projection='3d')
    CoorDataMask = np.arange(len(numLinks))[(np.sum((CoorData - CoorData[imagedHubs[0],:])**2,1)<((3*dthres)**2))]
    size = 60
    for i,(xs,ys,zs) in enumerate(CoorData[:,:]):
        if i in CoorDataMask:
            ax.scatter(xs, ys, zs, c='b', marker='o')
        if i in hubsList:
            ax.scatter(xs, ys, zs, c='r', marker='o',s=size)

if isInteractive:
    plt.show()

# end
