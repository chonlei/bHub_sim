# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 23:24:03 2017

@author: Chon Lei
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    speices = int(sys.argv[1])
    morphology = int(sys.argv[2])
except Exception:
    raise Exception("Please enter speices ID and morphology ID")
try:
    dcentre = int(sys.argv[3])
except Exception:
    dcentre = 40

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


# plot the islet
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
centreCoorData = np.mean(CoorData,0)
CoorDataMask = (np.sum((CoorData - centreCoorData)**2,1)<(dcentre**2))
for xs,ys,zs in CoorData[CoorDataMask,:]: #np.array(len(CoorData))[CoorDataMask]
    ax.scatter(xs, ys, zs, c='b', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')



# turn cm to usable coupling matrix in NEURON script
cm = np.tril(cm,-1)
#np.savetxt("CouplingMatrix-%s-%d.dat"%(coorDataFile,int(i*10),cm))

plt.show()