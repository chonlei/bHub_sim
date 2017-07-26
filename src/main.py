## Beta cell hub simulation
## Created by Chon Lei
## Last updated: August 2017

try:
    from neuron import h
except Exception:
    raise Exception("Please properly install NEURON package: http://www.neuron.yale.edu/neuron/download")
import numpy as np
import matplotlib.pylab as plt
import random as random
import sys
import modelSetup

## Define model and setup
model = 1
if model == 1:
    pathToModel = "../models/betacell_hermann2007_vFinal/"
    h('nrn_load_dll("%sx86_64/.libs/libnrnmech.so")'%pathToModel)
    h.load_file (pathToModel+"betacell.hoc")
    h.load_file (pathToModel+"betahub.hoc")
    h.load_file (pathToModel+"gapjunction.hoc")
    loadHetero = modelSetup.loadHeteroHermann2007
    setHetero = modelSetup.setHeteroHermann2007
    HetDict = modelSetup.HetDictHermann2007
## 
## Simulation of the electrophysiology of the beta-cell.
## Using package NEURON yale.
## 
## Created by Chon Lei
## The current version is based on the model from 
## M. Meyer-Hermann 2007 Biophysical Journal Vol. 93
## Last updated: 18/02/2017
## 
morphology = 1
if morphology == 1:
    CoupledMatrixFile = '../morphologies/mice/CouplingMatrix-mouse40-3-175.dat'
pyseed = 1
mode = 0
random.seed(pyseed)


## Import system setup files (.hoc files and system matrix)
CoupledMatrix = np.loadtxt(CoupledMatrixFile,delimiter=' ')
ncells = CoupledMatrix.shape[0]
if ncells != CoupledMatrix.shape[1]:
    raise Exception("CoupledMatrix invalid dimensions.")
try:
    pass
except AttributeError:
    raise Exception("Please make sure files has been compiled using \n$ nrnivmodl\n")
Total = (ncells*ncells)/2 - ncells # maximum number of gapjunctions that could be formed


## System set-up
print("*************************")
print("Starting system set-up...")

# Define beta hubs cells
numHubs = int(0.05*ncells)
temp = range(ncells)
random.shuffle(temp)
hubsList = temp[0:numHubs]

# Randomly pick some non-hubs cells
numNonHubsToPick = numHubs
temp = [i for i in range(ncells) if i not in hubsList]
random.shuffle(temp)
nonHubsToPickList = temp[0:numNonHubsToPick]

# Declare heterogeneity matrix
HetMatrix = np.zeros((len(HetDict)+1,ncells))
# Use a previously generated heterogeneity matrix
# Comment this out and use setHeteroHermann2007() to run another simulation
#HetMatrix = np.loadtxt('HetMatrix-mouse40-3.txt')

print "Defining cells..."
# Define as beta hub cell if in the hubsList
# Introduce heterogeneity to each defined cell
cell = []
iclamp_hubs = []
for i in range(ncells):
    if i not in hubsList:
        cell.append(h.betacell())
        setHetero(cell[i],HetMatrix,i)
        if (i in nonHubsToPickList) and (mode==2):
            iclamp_hubs.append(h.IClamp (0.5, sec = cell[i].soma) )
            iclamp_hubs[-1].delay = 0
            iclamp_hubs[-1].dur = 120000
            iclamp_hubs[-1].amp = -0.002
    else:
        cell.append(h.betahub())
        #loadHeteroHermann2007(cell[i],HetMatrix,i)
        cell[i].soma(0.5).nkatp_katp = -5.8
        if mode==1:
            # I clamp hubs to silence them, compare results from Johnston et al., 2016
            iclamp_hubs.append(h.IClamp (0.5, sec = cell[i].soma) )
            iclamp_hubs[-1].delay = 0
            iclamp_hubs[-1].dur = 120000
            #  all spiking: -0.0005; all stay -120mV: -0.005; all stay -72mV: -0.001; all stay -90mV: -0.002;
            iclamp_hubs[-1].amp = -0.002


print "Defining gap junction connections..."
gap = []
# 170pS (Zhang et al. 2008) Why the factor 1e-4
# Coupling too strong with 0.00017e-1 (potentially different units)
# empirically tested, using 3~5% of 0.00017e-1, weakly connected
for i in range(ncells):
    for j in range(ncells):
        if CoupledMatrix[i,j] > 0 and ((i in hubsList) or (j in hubsList)):
            #print i,j
            # strongly connect hubs with other cells, for testing purpose.
            gap.append(h.gapjunction(cell[i], cell[j], 0.5, 0.5, 0.5*0.00017e-1*CoupledMatrix[i,j]))
        elif CoupledMatrix[i,j] > 0: #and ((i not in nonHubsToPickList) or (j not in nonHubsToPickList)):
            gap.append(h.gapjunction(cell[i], cell[j], 0.5, 0.5, 0.05*0.00017e-1*CoupledMatrix[i,j]))


## External stimulation
"""
stimulus = h.IClamp (0.5, sec = cell[0].soma)
stimulus.delay = 100
stimulus.dur = 300
stimulus.amp = 0.5
"""


print "Creating recorder vectors..."
## System recorder initalisation
t = h.Vector()
t.record(h._ref_t)
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
h.tstop = 1e3
h.dt = 0.1
h.steps_per_ms = 1./h.dt
print("Starting main simulation...")
h.run()
print("Simulation completed! :)")
print("*************************")


## ******************************************** ##
## Exporting


## Visualisation
toVisual = False
if toVisual:
    t = np.array(t)
    membranepotential = np.array(vrec[0])
    plt.plot(t, membranepotential)
    
    membranepotential2 = np.array(vrec[102])
    plt.plot(t, membranepotential2)
    
    plt.show()
