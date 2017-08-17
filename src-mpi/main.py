#!/usr/bin/env python

import os
import sys
import simulator

from mpi4py import MPI


## description of current simulation
logmsg = "First testing of running multiple simulations using mpi+neuron\n\nUsing cubic lattice; common seed; not varying gGJ and gamma;silence one fixed hub in the image plane\n\nAdd hub one by one."


######
## Define model and setup
######
"""
model
gjmodel
morphology
species # 0: mouse; 1: human; 2: cubic lattice
pyseed 
isImitateExp # if True, simulate whole islet but only analyse imaged cells
mode # 0: WT; 1: silent hubs; 2: silent non hubs
silenceStart  # I clamp hubs to silence them, compare results from Johnston et al., 2016
silenceDur
silenceAmp #-100#mV  #-0.005#uA
pHubs  # percentage/fraction of hubs in islet (if <1) else number of hubs in islet (i.e. >1)
methodToPickHubs  # 0: random; 1: top GJ links; 2: bottom GJ links
whichHub # indix of imaged hub/non-hub to silence
ggap 
ggaphub 
pggaphubstd  # fraction of ggaphub as std
pggapstd  # fraction of ggap as std
gjtau 
dthres  # spatial cutoff distance to def GJ connection
isletsize # islet size of interest (None for whole islet)
hetVar  # it is % of mean's standard deviation
tstop # usually in [ms]
dt   # usually in [ms]
downSampling  # down sample the output -> output_timestep = dt*downSampling
tbatch  # split simulation into batches; same unit as tstop
"""
modelParam = {'model' : 2, \
              'gjmodel' : 1, \
              'morphology' : 1, \
              'species' : 2, \
              'pyseed' : 3, \
              'isImitateExp' : 1, \
              'mode' : 1, \
              'silenceStart' : 75e3, \
              'silenceDur' : 100e3, \
              'silenceAmp' : -0.005, \
              'pHubs' : 0.01, \
              'methodToPickHubs' : 0 , \
              'whichHub' : 0 , \
              'ggap' : 1.75/3.*1/6.*5.1*0.385*1e-4, \
              'ggaphub' : 0.8/3.*1/6.*5.1*0.385*1e-4, \
              'pggaphubstd' : 0, \
              'pggapstd' : 0.35, \
              'gjtau' : 100.0, \
              'dthres' : 17.5, \
              'isletsize' : 40 , \
              'hetVar' : 0.1, \
              'tstop' : 275e3, \
              'dt' : 0.1 , \
              'downSampling' : 1000, \
              'tbatch' : 5e3}

# model key word arguments
# model 1 default: {'beta':{} , 'betahub':{'hubkatp':-5.8}}
# model 2 default: {'beta':{'gkatp':(6.5,0.0) , 'useDistribution':None} , 'betahub':{'hubgkatp':10}}
modelParam['model_kwargs'] = { 'beta':{'gkatp':(4.0,7.0) , 'useDistribution':'sq'} , \
                               'betahub':{'hubgkatp':10} }

# setup output directory
#outputdir = modelSetup.outputMake()
try:
    outputdir = sys.argv[1]
    if not os.path.isdir(outputdir):
        print("Cannot set %s as output directory"%outputdir)
        outputdir = '../output/sim_temp/'
except Exception:
    print("Cannot set %s as output directory"%outputdir)
    outputdir = '../output/sim_temp/'
print("Output directory set as %s"%outputdir)
modelParam['parentout'] = outputdir


def main(comm,modelParam,outputdir):
    # Get rank and size
    size = comm.Get_size()
    rank = comm.Get_rank()
    modelParam['subidx'] = rank
    # setup what each sub simulation does
    modelParam['pHubs'] = 3*rank+3
    # check they are doing right thing
    #print("RANK %d of SIZE %d is doing the right job..."%(rank,size))
    simulator.main(modelParam)


if __name__ == "__main__":
    # log down some short note of current simulations
    # e.g. what is this simulation trying to test or what is changing.
    with open(os.path.join(outputdir,'logmsg.log'), 'w') as f:
        f.write("# SIMULATION DESCRIPTION\n\n")
        f.write(logmsg)
    # main simulation
    main(MPI.COMM_WORLD,modelParam,outputdir)

##
