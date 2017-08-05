# -*- coding: utf-8 -*-
"""
## 
## Simulation of the electrophysiology of the beta-cell and hubs effects.
## Using package NEURON yale.
##
## Try to split up the simulation into batches and compare it with one long simulation.
## This should be tested with 3 coupled cells.
## 
## Created by Chon Lei
## Last updated: August 2017
## 
"""

try:
    from neuron import h
except Exception:
    raise Exception("Please properly install NEURON package: http://www.neuron.yale.edu/neuron/download")
import numpy as np
import matplotlib.pylab as plt
import random as random
import os.path as path
import modelSetup


######
## Define model and setup
######
model = 2
gjmodel = 2
morphology = 3
ncells = 3
pyseed = 1
ggap = 1/6.*5.1*0.385*1e-4#0.5*0.00017e-1
ggaphub = 1/3.*1/6.*5.1*0.385*1e-4#1.0*0.00017e-1
gjtau = 160.0
hetVar = 0.05
tstop = 30e3  # usually in [ms]
dt = 0.1  # usually in [ms]
downSampling = 100  # down sample the output -> output_timestep = dt*downSampling
tbatch = 10e3 # split simulation into batches; same unit as tstop

# Create output directories and log file
outputidx, outputdir = modelSetup.outputSetup(model,morphology,pyseed,0,isTest=True)
outlog = path.join(outputdir, outputidx+'.log')
outCa = path.join(outputdir, 'Ca_'+outputidx)
outVm = path.join(outputdir, 'Vm_'+outputidx)
with open(outlog, 'w') as f:
    f.write('#model = %d \n#gjmodel = %d \n#morphology = %d \n#ncells = %d \n#pyseed = %d \n#ggap = %f \n#ggaphub = %f \n#gjtau = %f '%(model,gjmodel,morphology,ncells,pyseed,ggap,ggaphub,gjtau)+' \n')
    f.write('#hetVar = %f \n#tstop = %f \n#dt = %f \n#downSampling = %d \n#tbatch = %f \n\n'%(hetVar,tstop,dt,downSampling,tbatch))

if model == 1:
    ## Created by Chon Lei
    ## The current version is based on the model from 
    ## M. Meyer-Hermann 2007 Biophysical Journal Vol. 93
    ## Last updated: 18/02/2017
    pathToModel = "../models/betacell_hermann2007_vFinal/"
    modelSetup.setupHetDict(varp=hetVar)
    loadHetero = modelSetup.loadHeteroHermann2007
    setHetero = modelSetup.setHeteroHermann2007
    HetDict = modelSetup.HetDictHermann2007
    #h.load_file (pathToModel+"betahub.hoc")
    def defineBeta(cellList,i):
        # define beta cell
        cellList.append(h.betacell())
        setHetero(cellList[i],HetMatrix,i)
    def defineBetaHub(cellList,i):
        # define beta hub cell
        cellList.append(h.betahub())
        cellList[i].soma(0.5).nkatp_katp = -5.8
elif model == 2:
    ## Created by Chon Lei
    ## The current version is based on the model from
    ## M. Meyer-Hermann 2007 Biophysical Journal Vol. 93
    ## Last updated: 18/02/2017
    pathToModel = "../models/betacell_hermann2007_vMetabolic/"
    modelSetup.setupHetDict(varp=hetVar)
    loadHetero = modelSetup.loadHeteroHermann2007
    setHetero = modelSetup.setHeteroHermann2007
    HetDict = modelSetup.HetDictHermann2007
    def defineBeta(cellList,i):
        # define beta cell
        cellList.append(h.betacell())
        setHetero(cellList[i],HetMatrix,i)
        cellList[i].soma(0.5).gammatoset_katp = 7.0
    def defineBetaHub(cellList,i):
        # define beta hub cell
        cellList.append(h.betacell())
        cellList[i].soma(0.5).gammatoset_katp = 10.0

if gjmodel==1:
    pathToGJModel = "../models/gapjunction_pedersen2015/"
elif gjmodel==2:
    pathToGJModel = "../models/gapjunction_hermann2010/"

random.seed(pyseed)
np.random.seed(pyseed)
modelSetup.SetRandomSeed(pyseed)


######
## Import system setup files (.hoc files and system matrix)
######
CoupledMatrix = np.tril(np.ones((ncells,ncells)),-1) # just connect all cells
ncells = CoupledMatrix.shape[0]
nSpatialLinks = np.sum(CoupledMatrix+CoupledMatrix.T,1)
if ncells != CoupledMatrix.shape[1]:
    raise Exception("CoupledMatrix invalid dimensions.")
Total = (ncells*ncells)/2 - ncells # maximum number of gapjunctions that could be formed

try:
    # assumed using x64 bits; change if needed
    h('nrn_load_dll("%sx86_64/.libs/libnrnmech.so")'%pathToModel)
    h('nrn_load_dll("%sx86_64/.libs/libnrnmech.so")'%pathToGJModel)
    h.load_file (pathToModel+"betahub.hoc")
    h.load_file (pathToModel+"betacell.hoc")
    h.load_file (pathToGJModel+"gapjunction.hoc")
except Exception:
    raise Exception("Please make sure files has been compiled using \n$ nrnivmodl\n")

######
## System set-up
######
print("*************************")
print("Starting system set-up...")

# Declare heterogeneity matrix
HetMatrix = np.zeros((len(HetDict)+1,ncells))
# Use a previously generated heterogeneity matrix
# Comment this out and use setHeteroHermann2007() to run another simulation

print "Defining cells..."
# Define as beta hub cell if in the hubsList
# Introduce heterogeneity to each defined cell
cell = []
iclamp_hubs = []
for i in range(ncells):
    defineBeta(cell,i)
print HetMatrix

print "Defining gap junction connections..."
gap = []
# 170pS (Zhang et al. 2008) Why the factor 1e-4
# Coupling too strong with 0.00017e-1 (potentially different units)
# empirically tested, using 3~5% of 0.00017e-1, weakly connected
for i in range(ncells):
    for j in range(ncells):
        if CoupledMatrix[i,j] > 0:
            gap.append(h.gapjunction(cell[i], cell[j], 0.5, 0.5, ggap*CoupledMatrix[i,j],gjtau))

######
## External stimulation
######
"""
stimulus = h.IClamp (0.5, sec = cell[0].soma)
stimulus.delay = 100
stimulus.dur = 300
stimulus.amp = 0.5
"""

##TODO: do from here...


######
## System recorder initalisation
######
print("Creating recorder vectors...")
t = h.Vector()
t.record(h._ref_t)
vrec = []
carec = []
for i in range(ncells):
    vrec.append(h.Vector())
    carec.append(h.Vector())
    vrec[i].record(cell[i].soma(0.5)._ref_v)
    carec[i].record(cell[i].soma(0.5)._ref_cai)


######
## Main simulation
######
h.load_file("stdrun.hoc")
h.init()
# h.v_init = purkinjecelly.undershootpotential
#h.v_init = -70
h.tstop = tstop
h.dt = dt
h.steps_per_ms = 1./h.dt
print("Running main simulation...")
h.run()
print("Simulation completed! :)")
print("*************************")


## ******************************************** ##
######
## Exporting
######
print("Converting results for export...")
# Down sampling
carec = modelSetup.convertSimOutput(carec,downSampling)
vrec = modelSetup.convertSimOutput(vrec,downSampling)

# Output files
modelSetup.savedat(outCa,carec,'Ca',outlog,idx=None)
modelSetup.savedat(outVm,vrec,'Vm',outlog,idx=None)


"""
##TODO: check with one simulation; check duplication in successive saving.
temptstop = 0
nbatch = tstop%tbatch
for i in range(nbatch):
    h.frecord_init()
    h.continuerun(tstart)
"""

######
## Visualisation
######
toVisual = False
if toVisual:
    t = np.array(t)
    membranepotential = np.array(vrec[0])
    plt.plot(t, membranepotential)
    
    membranepotential2 = np.array(vrec[102])
    plt.plot(t, membranepotential2)
    
    plt.show()

##########
## End  ##
##########
