#!/usr/bin/env python

import os
import sys
import simulator

from mpi4py import MPI


## description of current simulation
logmsg = "First testing of running multiple simulations using mpi+neuron\n\nUsing mouse 1; (real) common seed; 1.5% hubs; silence n hubs in the islet\n\nn = 1, 2, 3...."#add each time 3 hubs from 12 to 3*3*16+12 hubs; silence first 10 hubs."#; not varying gGJ and gamma; 5% hubs; silence n hubs in the islet\n\nn = 1, 2, 3...."


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
ggap  # model 1,2: ~1/6.*5.1*0.385*1e-4; model 3: ~0.12 [nS]
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
try:
    pyseed = int(sys.argv[2])
except Exception:
    pyseed = 4
try:
    modee = int(sys.argv[3])
except Exception:
    modee = 1
try:
    archi = int(sys.argv[4])
except Exception:
    archi = 1
try:
    gjm = int(sys.argv[5])
except Exception:
    gjm = 1
try:
    gjhub = float(sys.argv[6])
except Exception:
    gjhub = 0.02
try:
    gjnonhub = float(sys.argv[7])
except Exception:
    gjnonhub = 0.02
try:
    mor = int(sys.argv[8])
except Exception:
    mor = 0
try:
    isAct = int(sys.argv[9])
except Exception:
    isAct = 1
try:
    stepsizer = int(sys.argv[10])
except Exception:
    stepsizer = 5
try:
    stepstart = int(sys.argv[11])
except Exception:
    stepstart = 0
try:
    base_glucose = float(sys.argv[12])
except Exception:
    base_glucose = 0
try:
    pgjhubstd = float(sys.argv[13])
except Exception:
    pgjhubstd = 0.2
try:
    pgjnonhubstd = float(sys.argv[14])
except Exception:
    pgjnonhubstd = 0.2
try:
    deltaPCaER = float(sys.argv[15])
    print('setting deltaPCaER to %g' % deltaPCaER)
except Exception:
    raise ValueError('Cannot read deltaPCaER')

modelParam = {'model' : 5, \
              'gjmodel' : gjm, \
              'morphology' : archi, \
              'species' : mor, \
              'pyseed' : pyseed, \
              'isImitateExp' : 1, \
              'mode' : modee, \
              'silenceStart' : 150e3, \
              'silenceDur' : 250e3, \
              'silenceAmp' : -100.0, \
              'pHubs' : 0.1, \
              'methodToPickHubs' : 0 , \
              'whichHub' : 0 , \
              'ggap' : gjnonhub, \
              'ggaphub' : gjhub, \
              'pggaphubstd' : pgjhubstd, \
              'pggapstd' : pgjnonhubstd, \
              'gjtau' : 400.0, \
              'p_connect': 1, \
              'dthres' : 17.5, \
              'isletsize' : 40 , \
              'hetVar' : 0.2, \
              'tstop' : 800e3, \
              'dt' : 0.1 , \
              'downSampling' : 1000, \
              'tbatch' : 5e3, \
              'deltaPCaER' : deltaPCaER}

# model key word arguments
# model 1 default: {'beta':{} , 'betahub':{'hubkatp':-5.8}}
# model 2 default: {'beta':{'gkatp':(6.5,0.0) , 'useDistribution':None} , 'betahub':{'hubgkatp':10}}
# model 3 default: {'beta':{'gkatp':(6.5,0.0) , 'useDistribution':None , 'applytime':5e3} , 'betahub':{'hubgkatp':10 , 'applytime':5e3}}
#modelParam['model_kwargs'] = {'beta':{'gkatp':(6.0,7.0) , 'useDistribution':'sq'} , 'betahub':{'hubgkatp':11.0}}
# modelParam['model_kwargs'] = {'beta':{'glu':(6.0,7.0) , 'useDistribution':'sq' , 'applytime':50e3} , \
#                               'betahub':{'hubglu':11.0 , 'applytime':50e3}}
modelParam['model_kwargs'] = {'beta':{'glu':(base_glucose,base_glucose+1.0) , 'useDistribution':'sq' , 'applytime':50e3} , \
                              'betahub':{'hubglu':11.0 , 'applytime':50e3}}

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
    if isAct == 1:
        modelParam['pHubs'] = 10*rank
        modelParam['silenceStart'] = 350e3
        modelParam['tstop'] = 360e3
    tempParam = stepstart + rank*stepsizer #rank+1 #(12*rank+30)/2 #
    # check they are doing right thing
    #print("RANK %d of SIZE %d is doing the right job..."%(rank,size))
    simulator.main(modelParam,tempParam=tempParam)


if __name__ == "__main__":
    # log down some short note of current simulations
    # e.g. what is this simulation trying to test or what is changing.
    with open(os.path.join(outputdir,'logmsg.log'), 'w') as f:
        f.write("# SIMULATION DESCRIPTION\n\n")
        f.write(logmsg)
    # main simulation
    main(MPI.COMM_WORLD,modelParam,outputdir)

##
