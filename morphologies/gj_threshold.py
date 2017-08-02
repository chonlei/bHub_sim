# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 23:24:03 2017

@author: Chon Lei
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings

try:
    speices = int(sys.argv[1])
    morphology = int(sys.argv[2])
except Exception:
    raise Exception("Please enter speices ID and morphology ID")
try:
    dcentre = np.float(sys.argv[3])
    isImitateExp = int(sys.argv[4])
except Exception:
    dcentre = 40
    isImitateExp = 0

# Define data file and parameters for testing
if speices==0:
    pathToFile = "./mouse/"
    coorDataFile = "Mouse 40-%d"%morphology
elif speices==1:
    #get human islet data
    pass
elif speices==2:
    pathToFile = "./cubic_lattice/"
    coorDataFile = "cubic%d"%morphology


cutoffTest = [5,7.5,10,12.5,15,17.5,20]
CoorData = np.loadtxt(pathToFile+coorDataFile+".txt")


def getImagedCellIdx(coorData, topDir=2, imageDepth=10, Ncells=100, method=0):
    """Get Experimentally trackable cell indices
    
    The function returns the cell indices that are experimentally trackable, 
    and there are two realisations of this, the first one simply takes the 
    cells that are within the imageDepth in the topDir direction; the second 
    method takes the cells that are at imageDepth distance underneath the cell 
    surface. Finally, you can also fix to pick the top Ncells number of cells. 
    
    Args:
        coorData (arr): The coordinate of the cells.
        topDir (int): The image direction/the direction to search cell. 
                      Pick: 0='x',1='y',2='z'. Default: 2.
        imageDepth (float): The distance of which can be imaged.
        Ncells (int): The number of cells to pick if using method=0. 
        method (int): The method to select cells.
                      0 = fixed number of cells to pick. (Default)
                      1 = imageDepth as distance from the top cell in the given
                          direction.
                      2 = cells that are at imageDepth distance underneath the 
                          islet surface.
    
    Return:
        cellIdx (arr): A list of cell indices that can be imaged.
    
    """
    cellIdx = []
    if method==0:
        # fixed number of cells to pick
        cellIdx = np.argsort(coorData[:,topDir])
        cellIdx = cellIdx[:Ncells] if Ncells>0 else cellIdx[Ncells:]
    elif method==1:
        # imageDepth as distance from the top cell in the given direction
        topCellCoor = np.max(coorData[:,topDir])
        isWithinDepth = coorData[:,topDir] > (topCellCoor-imageDepth)
        cellIdx = np.arange(len(coorData[:,topDir]))[isWithinDepth]
        if len(cellIdx)>Ncells:
            warnings.warn("Warning: returned cells are more than Ncells")
        elif len(cellIdx)<0.5*Ncells:
            warnings.warn("Warning: returned cells are less than half of Ncells")
    elif method==2:
        # cells that are at imageDepth distance underneath the islet surface
        # TBU 
        pass
    return cellIdx


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
cellImaged = getImagedCellIdx(CoorData)
if False:
    plt.figure(1)
    for i in cutoffTest:
        cm = genCoupleMatrix(CoorData,i)
        numLinks = np.sum(cm,1)
        countNumLinks = np.histogram(numLinks, bins=np.arange(numLinks.min(), numLinks.max()+1)-0.5)[0]
        plt.plot(np.arange(numLinks.min(), numLinks.max()), countNumLinks, '-s',label="$R_{th} = %g$"%i)
    
    plt.legend()
    plt.title("GJ spatial network distribution - %s"%coorDataFile)
    plt.xlabel("number of GJ to neighbouring cells")
    plt.ylabel("number of cells")
    #plt.savefig("SN_dist-%s.png"%coorDataFile)
    
    # turn cm to usable coupling matrix in NEURON script
    #cm = np.tril(cm,-1)
    #np.savetxt("CouplingMatrix-%s-%d.dat"%(coorDataFile,int(i*10),cm))
if True:
    cm = genCoupleMatrix(CoorData,17.5)
    plt.figure(3)
    cm1 = cm[cellImaged,:]
    numLinks = np.sum(cm1,1)
    countNumLinks = np.histogram(numLinks, bins=np.arange(numLinks.min(), numLinks.max()+1)-0.5)[0]
    plt.plot(np.arange(numLinks.min(), numLinks.max()), countNumLinks, '-s',label="all links")
    cm2 = cm[cellImaged,:][:,cellImaged]
    numLinks = np.sum(cm2,1)
    countNumLinks = np.histogram(numLinks, bins=np.arange(numLinks.min(), numLinks.max()+1)-0.5)[0]
    plt.plot(np.arange(numLinks.min(), numLinks.max()), countNumLinks, '-s',label="links amount imaged cells")
    plt.legend()
    plt.title("imaged cells GJ spatial network distribution - %s"%coorDataFile)
    plt.xlabel("number of GJ to neighbouring cells")
    plt.ylabel("number of cells")

# plot the islet
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
centreCoorData = np.mean(CoorData,0)
CoorDataMask = (np.sum((CoorData - centreCoorData)**2,1)<(dcentre**2))
print "imagedCells: ", cellImaged
for i,(xs,ys,zs) in enumerate(CoorData[CoorDataMask,:]): #np.array(len(CoorData))[CoorDataMask]
    ax.scatter(xs, ys, zs, c='b', marker='o')
    if i in cellImaged:
        ax.scatter(xs, ys, zs, c='r', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


plt.show()
