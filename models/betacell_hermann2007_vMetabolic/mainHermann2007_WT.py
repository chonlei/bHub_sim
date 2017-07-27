#! /usr/bin/python2.7
from neuron import h
import numpy as np
import matplotlib.pylab as plt
import random

## 
## Simulation of the electrophysiology of the beta-cell.
## Using package NEURON yale.
## 
## Created by Chon Lei
## The currently version is based on the model from 
## M. Meyer-Hermann 2007 Biophysical Journal Vol. 93
## Last updated: 04/04/2017
## 
## Included heterogeneity and hub cells properties 
## 


## Import system setup files (.hoc files and system matrix)
CoupledMatrix = np.loadtxt('CouplingMatrix-mouse40-3-175.dat',delimiter=' ')
ncells = CoupledMatrix.shape[0]
if ncells != CoupledMatrix.shape[1]:
    raise Exception("CoupledMatrix invalid dimensions.")

h.load_file ("betacell.hoc")
h.load_file ("betahub.hoc")
h.load_file("gapjunction.hoc")

Total = (ncells*ncells)/2 - ncells # maximum number of gapjunctions that could be formed



## ********************************************************** ##
## ********************  Some functions  ******************** ##
## ********************************************************** ##
#  to be sorted
#  to read, skip to the end of this session and continue

## Create a matrix to determine the parameters (e.g. conductance) of the cell-mechanism
#  For N-cells and M-parameters, define the matrix to be M-by-N matrix
#  
#  On the way, loop through the cell and assign to them.

HetIndexHermann2007 = {0:'rho_NaV', \
					1:'rho_NaK', \
					2:'rho_NCX', \
					3:'rho_KATP', \
					4:'rho_KV', \
					5:'rho_sKCa', \
					6:'rho_KCa', \
					7:'rho_CaL', \
					8:'rho_CaT', \
					9:'rho_PMCA'}

HetDictHermann2007 = {'rho_NaV':(1.4, 0.2*1.4), \
					'rho_NaK':(2000., 0.2*2000.), \
					'rho_NCX':(14., 0.2*14.), \
					'rho_KATP':(0.08, 0.2*0.08), \
					'rho_KV':(4.9, 0.2*4.9), \
					'rho_sKCa':(0.6, 0.2*0.6), \
					'rho_KCa':(0.4, 0.2*0.4), \
					'rho_CaL':(0.7, 0.2*0.7), \
					'rho_CaT':(0.15, 0.2*0.15), \
					'rho_PMCA':(1100., 0.2*1100.)}


def setHeteroHermann2007(b_cell,HetMatrix,i):
    # Passing b_cell, HetMatrix as a mutable object
    # b_cell: neuron h object constants function soma(0.5)
    # HetMatrix: matrix to store heterogeneity of cells
    # i: index (id) of b_cell
    b_cell.soma(0.5).rho_NaV_jna = HetMatrix[0,i] = np.random.normal(HetDictHermann2007['rho_NaV'][0],HetDictHermann2007['rho_NaV'][1])
    b_cell.soma(0.5).gbar_nap = HetMatrix[0,i]*b_cell.g_NaV_bar

    b_cell.soma(0.5).rho_NaK_jna = HetMatrix[1,i] = np.random.normal(HetDictHermann2007['rho_NaK'][0],HetDictHermann2007['rho_NaK'][1])
    b_cell.soma(0.5).rho_NaK_jk = HetMatrix[1,i]
    b_cell.soma(0.5).gbar_nak = HetMatrix[1,i]*b_cell.I_NaK_bar

    b_cell.soma(0.5).rho_NCX_jna = HetMatrix[2,i] = np.random.normal(HetDictHermann2007['rho_NCX'][0],HetDictHermann2007['rho_NCX'][1])
    b_cell.soma(0.5).rho_NCX_jca = HetMatrix[2,i]
    b_cell.soma(0.5).gbar_ncx = HetMatrix[2,i]*b_cell.I_NCX_bar

    b_cell.soma(0.5).rho_KATP_jk = HetMatrix[3,i] = np.random.normal(HetDictHermann2007['rho_KATP'][0],HetDictHermann2007['rho_KATP'][1])
    b_cell.soma(0.5).gbar_katp = HetMatrix[3,i]*b_cell.g_KATP_bar

    b_cell.soma(0.5).rho_KV_jk = HetMatrix[4,i] = np.random.normal(HetDictHermann2007['rho_KV'][0],HetDictHermann2007['rho_KV'][1])
    b_cell.soma(0.5).gbar_kv = HetMatrix[4,i]*b_cell.g_KV_bar

    b_cell.soma(0.5).rho_sKCa_jk = HetMatrix[5,i] = np.random.normal(HetDictHermann2007['rho_sKCa'][0],HetDictHermann2007['rho_sKCa'][1])
    b_cell.soma(0.5).gbar_skca = HetMatrix[5,i]*b_cell.g_sKCa_bar

    b_cell.soma(0.5).rho_KCa_jk = HetMatrix[6,i] = np.random.normal(HetDictHermann2007['rho_KCa'][0],HetDictHermann2007['rho_KCa'][1])
    b_cell.soma(0.5).gbar_kca = HetMatrix[6,i]*b_cell.g_KCa_bar

    b_cell.soma(0.5).rho_CaL_jca = HetMatrix[7,i] = np.random.normal(HetDictHermann2007['rho_CaL'][0],HetDictHermann2007['rho_CaL'][1])
    b_cell.soma(0.5).gbar_cal = HetMatrix[7,i]*b_cell.g_CaL_bar

    b_cell.soma(0.5).rho_CaT_jca = HetMatrix[8,i] = np.random.normal(HetDictHermann2007['rho_CaT'][0],HetDictHermann2007['rho_CaT'][1])
    b_cell.soma(0.5).gbar_cat = HetMatrix[8,i]*b_cell.g_CaT_bar

    b_cell.soma(0.5).rho_PMCA_jca = HetMatrix[9,i] = np.random.normal(HetDictHermann2007['rho_PMCA'][0],HetDictHermann2007['rho_PMCA'][1])
    b_cell.soma(0.5).gbar_pmca = HetMatrix[9,i]*b_cell.I_PMCA_bar

    b_cell.soma(0.5).nkatp_katp = HetMatrix[10,i] = np.random.normal(-6.5,0.3)
    b_cell.soma(0.5).kappa_KATP_jk = HetMatrix[10,i]


def loadHeteroHermann2007(b_cell,HetMatrix,i):
    # Passing b_cell, HetMatrix as a mutable object
    # b_cell: neuron h object constants function soma(0.5)
    # HetMatrix: matrix to assign heterogeneity of cells
    # i: index (id) of b_cell
    b_cell.soma(0.5).rho_NaV_jna = HetMatrix[0,i]
    b_cell.soma(0.5).gbar_nap = HetMatrix[0,i]*b_cell.g_NaV_bar

    b_cell.soma(0.5).rho_NaK_jna = HetMatrix[1,i]
    b_cell.soma(0.5).rho_NaK_jk = HetMatrix[1,i]
    b_cell.soma(0.5).gbar_nak = HetMatrix[1,i]*b_cell.I_NaK_bar

    b_cell.soma(0.5).rho_NCX_jna = HetMatrix[2,i]
    b_cell.soma(0.5).rho_NCX_jca = HetMatrix[2,i]
    b_cell.soma(0.5).gbar_ncx = HetMatrix[2,i]*b_cell.I_NCX_bar

    b_cell.soma(0.5).rho_KATP_jk = HetMatrix[3,i]
    b_cell.soma(0.5).gbar_katp = HetMatrix[3,i]*b_cell.g_KATP_bar

    b_cell.soma(0.5).rho_KV_jk = HetMatrix[4,i]
    b_cell.soma(0.5).gbar_kv = HetMatrix[4,i]*b_cell.g_KV_bar

    b_cell.soma(0.5).rho_sKCa_jk = HetMatrix[5,i]
    b_cell.soma(0.5).gbar_skca = HetMatrix[5,i]*b_cell.g_sKCa_bar

    b_cell.soma(0.5).rho_KCa_jk = HetMatrix[6,i]
    b_cell.soma(0.5).gbar_kca = HetMatrix[6,i]*b_cell.g_KCa_bar

    b_cell.soma(0.5).rho_CaL_jca = HetMatrix[7,i]
    b_cell.soma(0.5).gbar_cal = HetMatrix[7,i]*b_cell.g_CaL_bar

    b_cell.soma(0.5).rho_CaT_jca = HetMatrix[8,i]
    b_cell.soma(0.5).gbar_cat = HetMatrix[8,i]*b_cell.g_CaT_bar

    b_cell.soma(0.5).rho_PMCA_jca = HetMatrix[9,i]
    b_cell.soma(0.5).gbar_pmca = HetMatrix[9,i]*b_cell.I_PMCA_bar

    b_cell.soma(0.5).nkatp_katp = HetMatrix[10,i]
    b_cell.soma(0.5).kappa_KATP_jk = HetMatrix[10,i]


## Define the CoupledMatrix if provided a coordinate data
#  Procedure:
#  - load the coordinate
#  - define distance matrix
#    - by loop through the coordinates twice
#  - set cutoff distance cutoffConnect ~ 20 um
#  - convert it into binary matrix

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

## Tested
#A = np.array([[1, 2, 3],[3,2,1],[10,5,1],[100,1,2]])
#A = "testCoorData.txt"
#B = genCoupleMatrix(A)


## Define array to set b-cell hugs using cell index (id)
## Define b-cell hugs properties (potentially another .hoc)


## output CoupledMatrix, parameter matrix, hugs array
def setSimOutput(b_cells):
    # Return 3 lists: time, Vm, Ca_i
    # Except time, other outputs are list of variable of all cells
    # b_cells: a list containing b_cell objects
    n = len(b_cells)
    time = h.Vector()
    time.record(h._ref_t)
    vrec = []
    carec = []
    for i in range(n):
        vrec.append(h.Vector())
        carec.append(h.Vector())
        vrec[i].record(cell[i].soma(0.5)._ref_v)
        carec[i].record(cell[i].soma(0.5)._ref_cai)
    return time, vrec, carec

def convertSimOutput(listOutput,sparse=1):
    # Convert a list of h Vector object output to numpy array object
    # DO NOT USE if just want to convert a few of the cells in the whole list
    n = len(listOutput)
    for i in range(n):
        listOutput[i] = np.array(listOutput[i])[::sparse]
    return listOutput

## ********************************************************** ##
## ********************************************************** ##




## System set-up
print("*************************")
print("Starting system set-up...")

# Use a smaller number of cells to test
#ncells = 100
"""
array([  6.,   4.,   8.,  10.,   6.,   6.,  10.,  11.,  14.,  11.,  12.,
        14.,  10.,   5.,   7.,   9.,  12.,   8.,   8.,  11.,   8.,   9.,
        10.,  10.,  10.,   9.,   6.,   5.,   5.,   9.,   9.,  12.,  12.,
        11.,   3.,   7.,   4.,   5.,   8.,  11.,   9.,   6.,   5.,   4.,
         3.,   6.,  11.,   9.,   5.,  11.,   9.,  11.,  12.,  11.,  10.,
        10.,  11.,  10.,   7.,   8.,   5.,   5.,   3.,   3.,   8.,   7.,
         9.,   6.,   8.,   5.,   8.,   4.,   7.,   4.,   8.,   3.,   5.,
         9.,   8.,  10.,   9.,   8.,   6.,   4.,   8.,  11.,   8.,  11.,
         3.,   9.,  13.,  11.,  12.,  11.,   8.,  11.,  12.,  10.,   9.,
         6.])
"""

# Define beta hubs cells
numHubs = int(0.05*ncells)
temp = range(ncells)
random.shuffle(temp)
hubsList = temp[0:numHubs]
# Use a previously generated indices
# Comment this out to run another simulation
#hubsList = [64, 74, 34, 18, 11, 61, 98, 44, 94, 47] #[8] # k 14 # 
hubsList = [1025, 570, 1294, 81, 169, 659, 890, 1622, 1486, 1250, 247, 59, 595, 1526, 546, 1008, 1629, 1748, 923, 872, 742, 635, 920, 977, 1333, 867, 1438, 434, 524, 1053, 1235, 420, 718, 1042, 835, 399, 775, 1275, 1573, 148, 407, 365, 1240, 515, 523, 452, 391, 1743, 1049, 753, 1463, 388, 1502, 448, 166, 1718, 687, 1090, 1536, 254, 1219, 1490, 683, 0, 584, 1701, 1747, 11, 1424, 1350, 377, 1062, 1031, 1026, 1266, 1652, 367, 1181, 669, 550, 643, 1137, 1349, 1411, 141, 1247, 95, 1082] # same as strong GJ
print( "hubs list: ", hubsList )


# Randomly pick some non-hubs cells
numNonHubsToPick = numHubs
temp = [i for i in range(ncells) if i not in hubsList]
random.shuffle(temp)
nonHubsToPickList = temp[0:numNonHubsToPick]
# Use a previously generated indices
# Comment this out to run another simulation
nonHubsToPickList = [949, 484, 1372, 974, 1160, 804, 63, 442, 1303, 302, 192, 1726, 539, 313, 851, 203, 494, 552, 482, 1081, 1166, 1110, 697, 665, 1663, 764, 1758, 979, 535, 929, 1676, 904, 1001, 1340, 1503, 103, 1309, 1686, 402, 1492, 548, 847, 65, 1679, 1302, 28, 1626, 1678, 1390, 1405, 217, 1224, 845, 160, 36, 720, 496, 120, 739, 1633, 492, 1671, 1106, 1024, 1209, 156, 640, 1693, 269, 1508, 1242, 1119, 512, 1713, 1432, 922, 557, 1107, 676, 271, 1195, 409, 1725, 628, 1610, 1229, 1498, 341] #[28, 35, 89, 38, 27, 56, 65, 5, 6, 17]  # [15, 83, 60, 86, 28, 35, 58, 63, 46, 92] #
print( "nonhubs to pick: ", nonHubsToPickList )


# heterogeneity matrix
HetMatrix = np.zeros((len(HetDictHermann2007)+1,ncells))
# Use a previously generated heterogeneity matrix
# Comment this out and use setHeteroHermann2007() to run another simulation
HetMatrix = np.loadtxt('HetMatrix-mouse40-3.txt')

print( "Defining cells..." )
# Define as beta hub cell if in the hubsList
# Introduce heterogeneity to each defined cell
cell = []
iclamp_hubs = []
for i in range(ncells):
    if i not in hubsList:
        cell.append(h.betacell())
        cell[i].soma(0.5).gammatoset_katp = 7.0
        loadHeteroHermann2007(cell[i],HetMatrix,i)
        '''if i in nonHubsToPickList:
            iclamp_hubs.append(h.IClamp (0.5, sec = cell[i].soma) )
            iclamp_hubs[-1].delay = 0
            iclamp_hubs[-1].dur = 120000
            iclamp_hubs[-1].amp = -0.002'''
    else:
        cell.append(h.betacell())
        #loadHeteroHermann2007(cell[i],HetMatrix,i)
        cell[i].soma(0.5).gammatoset_katp = 10.0
        if i in hubsList[:len(hubsList)]:
        # I clamp hubs to silence them, compare results from Johnston et al., 2016
            iclamp_hubs.append(h.IClamp (0.5, sec = cell[i].soma) )
            iclamp_hubs[-1].delay = 0
            iclamp_hubs[-1].dur = 120000
        #  all spiking: -0.0005; all stay -120mV: -0.005; all stay -72mV: -0.001; all stay -90mV: -0.002; 
            iclamp_hubs[-1].amp = -0.002
        

print( "Defining gap junction connections..." )
gap = []
# 170pS (Zhang et al. 2008) Why the factor 1e-4
# Coupling too strong with 0.00017e-1 (potentially different units)
# empirically tested, using 3~5% of 0.00017e-1, weakly connected
for i in range(ncells):
    for j in range(ncells):
        if CoupledMatrix[i,j] > 0 and ((i in hubsList) or (j in hubsList)):
            #print( i,j )
            # strongly connect hubs with other cells, for testing purpose.
            gap.append(h.gapjunction(cell[i], cell[j], 0.5, 0.5, 1.0*0.00017e-1*CoupledMatrix[i,j]))
        elif CoupledMatrix[i,j] > 0: #and ((i not in nonHubsToPickList) or (j not in nonHubsToPickList)):
            gap.append(h.gapjunction(cell[i], cell[j], 0.5, 0.5, 0.5*0.00017e-1*CoupledMatrix[i,j]))




## External stimulation
"""
stimulus = h.IClamp (0.5, sec = cell[0].soma)
stimulus.delay = 100
stimulus.dur = 300
stimulus.amp = 0.5
"""

print( "Creating recorder vectors..." )
## System recorder initalisation
time = h.Vector()
time.record(h._ref_t)
vrec = []
carec = []
for i in range(ncells):
    vrec.append(h.Vector())
    carec.append(h.Vector())
    vrec[i].record(cell[i].soma(0.5)._ref_v)
    carec[i].record(cell[i].soma(0.5)._ref_cai)



## Main simulation
h.load_file("stdrun.hoc")
h.init()
# h.v_init = purkinjecelly.undershootpotential
#h.v_init = -70
h.tstop = 50e3
h.dt = 0.1
h.steps_per_ms = 1./h.dt
print("Starting main simulation...")
h.run()
print("Simulation completed! :)")
print("*************************")


## ******************************************** ##
## Export

savename = "p63_210417_1hub_"

np.savetxt('Het20pMatrix-mouse40-3.txt',HetMatrix)
carec = convertSimOutput(carec,100)
vrec = convertSimOutput(vrec,100)
#np.savetxt('Ca_mnonhubs_10hubs_90s.txt',carec)
#np.savetxt('Vm_mnonhubs_10hubs_90s.txt',vrec)

print( "Exporting Ca traces..." )

# Export using np.memmap
from tempfile import mkdtemp
import os.path as path
#filename = path.join(mkdtemp(), '%sca_dynamics.dat'%savename)
filename = path.join('../out/', 'Ca_WThet20p_5pGLhubs_50s_sGJ05_%dby%d.dat'%(len(carec),len(carec[0])))

fp = np.memmap(filename, dtype="float64", mode='w+', shape=(len(carec),len(carec[0])))
fp[:] = carec[:]
if fp.filename == path.abspath(filename):
	print( "shape: (", len(carec), ", ", len(carec[0]), ")" )
	print("Successfully exported Ca_dynamics to: %s"%filename)
	del fp
else:
	print("Cannot write to address: %s"%filename)
	del fp
	print("Writing to current path instead...")
	np.savetxt('carec.txt',carec)
	print("Successfully exported Ca_dynamics to current path.")
# newfp = np.memmap(filename, dtype='float64', mode='r', shape=(#cells,#timepoints)) # to read


print( "Exporting Vm traces..." )

filename = path.join('../out/', 'Vm_WThet20p_5pGLhubs_50s_sGJ05_%dby%d.dat'%(len(vrec),len(vrec[0])))

fp = np.memmap(filename, dtype="float64", mode='w+', shape=(len(vrec),len(vrec[0])))
fp[:] = vrec[:]
if fp.filename == path.abspath(filename):
	print( "shape: (", len(vrec), ", ", len(vrec[0]), ")" )
	print("Successfully exported Vm to: %s"%filename)
	del fp
else:
	print("Cannot write to address: %s"%filename)
	del fp
	print("Writing to current path instead...")
	np.savetxt('vrec.txt',vrec)
	print("Successfully exported Vm to current path.")



## ******************************************** ##


'''
## Visualisation
time = np.array(time)
membranepotential = np.array(vrec[hubsList[0]])
plt.plot(time, membranepotential-70)
for i in range(10):
	membranepotential = np.array(vrec[i])
	plt.plot(time, membranepotential+70*i)

plt.figure(2)
Caave = np.zeros(carec[0].shape)
for i in range(ncells):
	if (i not in hubsList):
		plt.plot(time, carec[i], 'k', alpha=0.5)
		Caave += carec[i]
plt.plot(time, carec[i], 'k', alpha=0.5, label="followers")
for i in hubsList:
	plt.plot(time, carec[i], 'r', linewidth=2)
	Caave += carec[i]
plt.plot(time, carec[i], 'r', linewidth=2, label="hubs")
plt.xlabel("t [ms]")
plt.ylabel("[Ca]i [mM]")
plt.legend()
#plt.savefig("%sca.png"%savename)

plt.figure(3)
Vave = np.zeros(vrec[0].shape)
for i in range(ncells):
	if (i not in hubsList):
		plt.plot(time, vrec[i], 'k', alpha=0.5)
		Vave += vrec[i]
plt.plot(time, vrec[i], 'k', alpha=0.5, label="followers")
for i in hubsList:
	plt.plot(time, vrec[i], 'r', linewidth=2)
	Vave += vrec[i]
plt.plot(time, vrec[i], 'r', linewidth=2, label="hubs")
plt.xlabel("t [ms]")
plt.ylabel("V [mV]")
plt.legend()
#plt.savefig("%sv.png"%savename)

Vave = Vave/ncells
Caave = Caave/ncells

plt.figure(4)
plt.plot(time, Vave)
plt.xlabel("t [ms]")
plt.ylabel("V [mV]")
plt.title("Averaged out all cells")

plt.figure(5)
plt.plot(time, Caave)
plt.xlabel("t [ms]")
plt.ylabel("[Ca]i [mM]")
plt.title("Averaged out all cells")

#np.savetxt('%sv_ave.txt'%savename,Vave)
#np.savetxt('%sca_ave.txt'%savename,Caave)

plt.show()
'''
