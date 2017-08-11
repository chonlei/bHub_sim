# -*- coding: utf-8 -*-
"""
## 
## Simulation of the electrophysiology of the beta-cell and hubs effects.
## Using package NEURON yale.
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
morphology = 1
species = 0  # 0: mouse; 1: human; 2: cubic lattice
pyseed = 2
isImitateExp = 1  # if True, simulate whole islet but only analyse imaged cells
mode = 1  # 0: WT; 1: silent hubs; 2: silent non hubs
silenceStart = 30e3#50e3
silenceDur = 250e3
silenceAmp = -0.005
pHubs = 0.01  # percentage/fraction of hubs in islet
##TODO need to do methodToPickHubs
methodToPickHubs = 0  # 0: random; 1: top GJ links; 2: bottom GJ links
ggap = 1/6.*5.1*0.385*1e-4#0.5*0.00017e-1
ggaphub = 1/3.*1/6.*5.1*0.385*1e-4#1.0*0.00017e-1
gjtau = 160.0
dthres = 17.5  # spatial cutoff distance to def GJ connection
isletsize = 40  # islet size of interest (None for whole islet)
hetVar = 0.05
tstop = 5e3#350e3  # usually in [ms]
dt = 0.1  # usually in [ms]; simulation dt
downSampling = 1  # down sample the output -> output_timestep = Dt*downSampling
tbatch = 1e3 # split simulation into batches; same unit as tstop
Dt = 100.0  # neuron recording vector down sample; neuron output dt

if isImitateExp:
    isletsize = None # force to use whole islet

# Create output directories and log file
outputidx, outputdir = modelSetup.outputSetup(model,morphology,pyseed,mode)
outlog = path.join(outputdir, outputidx+'.log')
outCa = path.join(outputdir, 'Ca_'+outputidx)
outVm = path.join(outputdir, 'Vm_'+outputidx)
with open(outlog, 'w') as f:
    f.write('#model = %d \n#gjmodel = %d \n#morphology = %d \n#species = %d \n#pyseed = %d \n#isImitateExp = %d \n#mode = %d \n#silenceStart = %f \n#silenceDur = %f \n#silenceAmp = %f \n#pHubs = %f \n#methodToPickHubs = %d \n#ggap = %f \n#ggaphub = %f \n#gjtau = %f \n#dthres = %f \n#isletsize = '%(model,gjmodel,morphology,species,pyseed,isImitateExp,mode,silenceStart,silenceDur,silenceAmp,pHubs,methodToPickHubs,ggap,ggaphub,gjtau,dthres)+str(isletsize)+' \n')
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

if species==0:
    #pathToCoupledMatrix = '../morphologies/mouse/CouplingMatrix-mouse40-3-175.dat'
    #pathToCoupledMatrix = '../morphologies/mouse/CouplingMatrixMouse403.dat'
    pathToMorphology = "../morphologies/mouse/Mouse 40-%d.txt"%morphology
elif species==1:
    pathToMorphology = ""
elif species==2:
    pathToMorphology = "../morphologies/cubic_lattice/cubic%d.txt"%morphology

random.seed(pyseed)
np.random.seed(pyseed)
modelSetup.SetRandomSeed(pyseed)


######
## Import system setup files (.hoc files and system matrix)
######
#CoupledMatrix = np.loadtxt(pathToCoupledMatrix,delimiter=' ')#[-100:,-100:] #only taking last 100 cells
CoorData = np.loadtxt(pathToMorphology)
# process CoorData to be acceptable format in modelSetup.genCoupleMatrix()
CoorData = CoorData[CoorData[:,0]==11][:,1:4]
CoupledMatrix = modelSetup.genCoupleMatrix(CoorData,dthres,isletsize,True)
ncells = CoupledMatrix.shape[0]
nSpatialLinks = np.sum(CoupledMatrix+CoupledMatrix.T,1)
if ncells != CoupledMatrix.shape[1]:
    raise Exception("CoupledMatrix invalid dimensions.")
Total = (ncells*ncells)/2 - ncells # maximum number of gapjunctions that could be formed

if isImitateExp==1:
    imagedCells = modelSetup.getImagedCellIdx(CoorData,topDir=2,imageDepth=10,Ncells=100,method=0)
with open(outlog, 'a') as f:
    f.write('#imagedCells = ')
    f.write(','.join(map(str, imagedCells)))
    f.write('\n\n')

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

# Define beta hubs cells
numHubs = int(pHubs*ncells)
temp = range(ncells)
random.shuffle(temp)
hubsList = temp[0:numHubs]
imagedHubs = list(set(hubsList).intersection(imagedCells))
#hubsList.remove(imagedHubs[1])
#imagedHubs.remove(imagedHubs[1])
if len(imagedHubs) < int(pHubs*len(imagedCells)):
    nMorehubs = int(pHubs*len(imagedCells)) - len(imagedHubs)
    tempCellsToPick = [x for x in imagedCells if x not in imagedHubs]
    random.shuffle(tempCellsToPick)
    hubsList += tempCellsToPick[0:nMorehubs]
    imagedHubs += tempCellsToPick[0:nMorehubs]
# Use a previously generated indices
# Uncomment this to set specific hubs
#hubsList = [self_def_list]
print(hubsList)
with open(outlog, 'a') as f:
    f.write('#hubsList = \n')
    f.write(','.join(map(str, hubsList)))
    f.write('\n\n')
    f.write('#imagedHubs = ')
    f.write(','.join(map(str, imagedHubs)))
    f.write('\n\n')

# Randomly pick some non-hubs cells
numNonHubsToPick = numHubs
temp = [i for i in range(ncells) if i not in hubsList]
random.shuffle(temp)
nonHubsToPickList = temp[0:numNonHubsToPick]
imagedNonHubs = []
if isImitateExp:
    temp = [i for i in imagedCells if i not in imagedHubs]
    random.shuffle(temp)
    imagedNonHubs = temp[0:int(pHubs*len(imagedCells))]
print(nonHubsToPickList)
with open(outlog, 'a') as f:
    f.write('#nonHubsToPickList = ')
    f.write(','.join(map(str, nonHubsToPickList)))
    f.write('\n\n')
    f.write('#imagedNonHubs = ')
    f.write(','.join(map(str, imagedNonHubs)))
    f.write('\n\n')

# Declare heterogeneity matrix
HetMatrix = np.zeros((len(HetDict)+1,ncells))
# Use a previously generated heterogeneity matrix
# Comment this out and use setHeteroHermann2007() to run another simulation
#HetMatrix = np.loadtxt('HetMatrix-mouse40-3.txt')

##TODO: Add function to silent cells
def silenceCell(iclampList,cell,delay=250e3,dur=250e3,amp=-0.005):
    # iclamp cell to imitate cell silencing in experiments
    # all spiking: -0.0005; all stay -120mV: -0.005; all stay -72mV: -0.001; all stay -90mV: -0.002;
    iclampList.append(h.IClamp (0.5, sec = cell.soma) )
    iclampList[-1].delay = delay
    iclampList[-1].dur = dur
    iclampList[-1].amp = amp
print "Defining cells..."
# Define as beta hub cell if in the hubsList
# Introduce heterogeneity to each defined cell
cell = []
iclamp_hubs = []
for i in range(ncells):
    if i not in hubsList:
        #cell.append(h.betacell())
        #setHetero(cell[i],HetMatrix,i)
        defineBeta(cell,i)
        if isImitateExp:
            if (i == imagedNonHubs[0]) and (mode==2):
                print "silencing cell ",i
                with open(outlog, 'a') as f:
                    f.write('#silencedCell = %d\n'%i)
                    f.write('#cell%d_nSpatialLinks = %d\n'%(i,nSpatialLinks[i]))
                silenceCell(iclamp_hubs,cell[i],silenceStart,silenceDur,silenceAmp)
        else:
            if (i in nonHubsToPickList) and (mode==2):
                print "silencing cell ",i
                with open(outlog, 'a') as f:
                    f.write('#silencedCell = %d\n'%i)
                    f.write('#cell%d_nSpatialLinks = %d\n'%(i,nSpatialLinks[i]))
                silenceCell(iclamp_hubs,cell[i],silenceStart,silenceDur,silenceAmp)
    else:
        #cell.append(h.betahub())
        ##loadHetero(cell[i],HetMatrix,i) # no heterogenity for betahubs
        #cell[i].soma(0.5).nkatp_katp = -5.8
        defineBetaHub(cell,i)
        if isImitateExp:
            if mode==1 and (i == imagedHubs[0]):
                # I clamp hubs to silence them, compare results from Johnston et al., 2016
                print "silencing cell ",i
                with open(outlog, 'a') as f:
                    f.write('#silencedCell = %d\n'%i)
                    f.write('#cell%d_nSpatialLinks = %d\n'%(i,nSpatialLinks[i]))
                silenceCell(iclamp_hubs,cell[i],silenceStart,silenceDur,silenceAmp)
        else:
            if mode==1:
                # I clamp hubs to silence them, compare results from Johnston et al., 2016
                print "silencing cell ",i
                with open(outlog, 'a') as f:
                    f.write('#silencedCell = %d\n'%i)
                    f.write('#cell%d_nSpatialLinks = %d\n'%(i,nSpatialLinks[i]))
                silenceCell(iclamp_hubs,cell[i],silenceStart,silenceDur,silenceAmp)


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
            gap.append(h.gapjunction(cell[i], cell[j], 0.5, 0.5, ggaphub*CoupledMatrix[i,j],gjtau))
        elif CoupledMatrix[i,j] > 0: #and ((i not in nonHubsToPickList) or (j not in nonHubsToPickList)):
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


######
## System recorder initalisation
######
print("Creating recorder vectors...")
t = h.Vector()
t.record(h._ref_t)
vrec = []
carec = []
for i in range(ncells):
    vrec.append(h.Vector(int(tbatch/Dt)))
    carec.append(h.Vector(int(tbatch/Dt)))
    vrec[i].record(cell[i].soma(0.5)._ref_v,Dt)
    carec[i].record(cell[i].soma(0.5)._ref_cai,Dt)


######
## Main simulation
######
h.load_file("stdrun.hoc")
h.init()
h.dt = dt
h.steps_per_ms = 1./h.dt
print("Running main simulation...")

temptstop = 0  # tstop for current batch
nbatch = int(tstop/tbatch)  # split simulation into nbatch
print("Dividing simulation into %d batches..."%nbatch)
tremain = tstop%tbatch  # remaining simulation time after nbatch
for i in xrange(nbatch):
    temptstop += tbatch  # tstop for current batch
    if False:
        for ii in range(ncells):
            vrec[ii].resize(0)
            carec[ii].resize(0)
            vrec[ii] = h.Vector(int(tbatch/Dt))
            carec[ii] = h.Vector(int(tbatch/Dt))
            vrec[ii].record(cell[ii].soma(0.5)._ref_v,Dt)
            carec[ii].record(cell[ii].soma(0.5)._ref_cai,Dt)
    h.frecord_init()  # reuse all recording vectors
    h.continuerun(temptstop)
    # exporting Ca time series
    tosave = modelSetup.convertSimOutput(carec,downSampling,reuse=True)
    modelSetup.savedat(outCa,tosave,'Ca',outlog,idx=i)
    # exporting Vm time series
    tosave = modelSetup.convertSimOutput(vrec,downSampling,reuse=True)
    modelSetup.savedat(outVm,tosave,'Vm',outlog,idx=i)
    print("Finished section %d out of %d."%(i+1,nbatch))
if tremain > 0:
    print("Running final section...")
    h.frecord_init()  # reuse all recording vectors
    h.continuerun(tstop)  # run until the end
    # exporting Ca time series
    tosave = modelSetup.convertSimOutput(carec,downSampling,reuse=True)
    modelSetup.savedat(outCa,tosave,'Ca',outlog,idx=i+1)
    # exporting Vm time series
    tosave = modelSetup.convertSimOutput(vrec,downSampling,reuse=True)
    modelSetup.savedat(outVm,tosave,'Vm',outlog,idx=i+1)
print("Simulation completed! :)")
print("*************************")


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
