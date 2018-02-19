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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import random as random
import os.path as path
import sys
sys.path.append("../src/")
import modelSetup

import pickle


######
## Define model and setup; for testing only
######
"""
model
gjmodel
morphology
species # 0: mouse; 1: human; 2: cubic lattice
pyseed 
isImitateExp # if True, simulate whole islet but only analyse imaged cells
mode # 0: WT; 1: silent hubs; 2: silent non hubs; 3: inject current to delta cells
silenceStart  # I clamp hubs to silence them, compare results from Johnston et al., 2016
silenceDur
silenceAmp #-100#mV  #-0.005#uA # if mode=3, this will apply to delta cells and use +0.09#uA
pHubs  # percentage/fraction of hubs in islet (if <1) else number of hubs in islet (i.e. >1)
methodToPickHubs  # 0: random; 1: top GJ links; 2: bottom GJ links; 3: next to delta cells
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
modelParam = {'model' : 5, \
              'gjmodel' : 1, \
              'morphology' : 5, \
              'species' : 1, \
              'pyseed' : 11, \
              'isImitateExp' : 1, \
              'mode' : 3, \
              'silenceStart' : 150e3, \
              'silenceDur' : 250e3, \
              'silenceAmp' : 0.09, \
              'pHubs' : 0.1, \
              'methodToPickHubs' : 3 , \
              'whichHub' : 0 , \
              'ggap' : 0.01, \
              'ggaphub' : 0.05, \
              'pggaphubstd' : 0.1, \
              'pggapstd' : 0.7, \
              'gjtau' : 400.0, \
              'p_connect': 1, \
              'dthres' : 17.5, \
              'isletsize' : 40 , \
              'hetVar' : 0.2, \
              'tstop' : 600e3, \
              'dt' : 0.1 , \
              'downSampling' : 1000, \
              'tbatch' : 5e3}

#modelParam['model_kwargs'] = {'beta':{'gamma':(0.985,0.0) , 'useDistribution':None , 'applytime':0e3} , 'betahub':{'hubgamma':1.0 , 'applytime':0e3}}
#modelParam['model_kwargs'] = {'beta':{'gkatp':(6.5,6.5) , 'useDistribution':'sq' , 'applytime':0e3} , 'betahub':{'hubgkatp':11.0 , 'applytime':0e3}}
modelParam['model_kwargs'] = {'beta':{'glu':(6.0,7.0) , 'useDistribution':'sq' , 'applytime':50e3} , \
                              'betahub':{'hubglu':11.0 , 'applytime':50e3}}


def main(modelParam=modelParam, hubsList_temp=[]):
    ######
    ## Define model and setup
    ######
    model = modelParam['model' ]
    gjmodel = modelParam['gjmodel']
    morphology = modelParam['morphology']
    species = modelParam['species']
    pyseed = modelParam['pyseed']
    isImitateExp = modelParam['isImitateExp']
    mode = modelParam['mode']
    silenceStart = modelParam['silenceStart']
    silenceDur = modelParam['silenceDur']
    silenceAmp = modelParam['silenceAmp']
    pHubs = modelParam['pHubs']
    ##TODO need to do methodToPickHubs
    methodToPickHubs = modelParam['methodToPickHubs']
    whichHub = modelParam['whichHub']
    ggap = modelParam['ggap']
    ggaphub = modelParam['ggaphub']
    try:
        pggaphubstd = modelParam['pggaphubstd']
    except Exception:
        pggaphubstd = 0.0
    try:
        pggapstd = modelParam['pggapstd']
    except Exception:
        pggapstd = 0.0
    gjtau = modelParam['gjtau']
    dthres = modelParam['dthres']
    isletsize = modelParam['isletsize']
    hetVar = modelParam['hetVar']
    tstop = modelParam['tstop']
    dt = modelParam['dt']
    downSampling = modelParam['downSampling']
    tbatch = modelParam['tbatch']
    try:
        # model 1 default: {'beta':{} , 'betahub':{'hubkatp':-5.8}}
        # model 2 default: {'beta':{'gkatp':(6.5,0.0) , 'useDistribution':None} , 'betahub':{'hubgkatp':10}}
        # model 3 default: {'beta':{'gkatp':(6.5,0.0) , 'useDistribution':None , 'applytime':5e3} , 'betahub':{'hubgkatp':10 , 'applytime':5e3}}
        # model 4 default: {'beta':{'gamma':(0.5,0.0) , 'useDistribution':None , 'applytime':5e3} , 'betahub':{'hubgamma':1.0 , 'applytime':5e3}}
        model_kwargs = modelParam['model_kwargs']
    except:
        model_kwargs = { 'beta':{} , 'betahub':{} }
    try:
        p_connect = modelParam['p_connect']
    except:
        p_connect = 1.0

    if isImitateExp:
        isletsize = None # force to use whole islet

    # Create output directories and log file
    try:
        parentout = modelParam['parentout']
        subidx = modelParam['subidx']
        outputidx, outputdir = modelSetup.outputSetup_sub(model,morphology,pyseed,mode,parentout,subidx)
    except Exception:
        outputidx, outputdir = modelSetup.outputSetup(model,morphology,pyseed,mode, isTest=True)
    outlog = path.join(outputdir, outputidx+'.log')
    outCa = path.join(outputdir, 'Ca_'+outputidx)
    outVm = path.join(outputdir, 'Vm_'+outputidx)
    outDeltaCa = path.join(outputdir, 'DeltaCa_'+outputidx)
    outDeltaVm = path.join(outputdir, 'DeltaVm_'+outputidx)
    # save modelParam
    with open(path.join(outputdir,'modelParam.pkl'), 'wb') as f:
        pickle.dump(modelParam, f, pickle.HIGHEST_PROTOCOL)
    # log model setup as txt file
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
        def defineBeta(cellList,i,**kwargs):
            # define beta cell
            cellList.append(h.betacell())
            if hetVar>0:
                setHetero(cellList[i],HetMatrix,i)
        def defineBetaHub(cellList,i,hubkatp=-5.8):
            # define beta hub cell
            cellList.append(h.betahub())
            cellList[i].soma(0.5).nkatp_katp = hubskatp
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
        def defineBeta(cellList,i,gkatp=(6.5,0.0),useDistribution=None):
            # define beta cell
            cellList.append(h.betacell())
            if hetVar>0:
                setHetero(cellList[i],HetMatrix,i)
            if useDistribution == None:
                cellList[i].soma(0.5).gammatoset_katp = gkatp[0]
            elif useDistribution == 'sq':
                cellList[i].soma(0.5).gammatoset_katp = np.random.uniform(gkatp[0],gkatp[1])
            elif useDistribution == 'normal':
                cellList[i].soma(0.5).gammatoset_katp = gkatp[0] + np.random.normal(0.0,1.0)*np.sqrt(gkatp[1])
            else:
                cellList[i].soma(0.5).gammatoset_katp = gkatp[0]
        def defineBetaHub(cellList,i,hubgkatp=10.0):
            # define beta hub cell
            cellList.append(h.betacell())
            cellList[i].soma(0.5).gammatoset_katp = hubgkatp
    elif model==3:
        ## Created by Chon Lei
        ## The current version is based on the model from
        ## Cha et al. 2011 The Journal of General Physiology
        ## Last updated: 18/08/2017
        pathToModel = "../models/betacell_cha2011_vMetabolic/"
        #TODO
        modelSetup.setupHetDict(varp=hetVar)
        #loadHetero = modelSetup.setHeteroCha2011
        setHetero = modelSetup.setHeteroCha2011
        setOrigin = modelSetup.setOriginCha2011
        HetDict = modelSetup.HetDictCha2011
        def defineBeta(cellList,i,gkatp=(6.0,0.0),useDistribution=None,applytime=5e3):
            # define beta cell
            cellList.append(h.betacell())
            cellList[i].soma(0.5).gammaapplytime_bcellcha = applytime
            if hetVar>0:
                setHetero(cellList[i],HetMatrix,i)
            else:
                setOrigin(cellList[i])
            if useDistribution == None:
                cellList[i].soma(0.5).gammatoset_bcellcha = gkatp[0]
            elif useDistribution == 'sq':
                cellList[i].soma(0.5).gammatoset_bcellcha = np.random.uniform(gkatp[0],gkatp[1])
            elif useDistribution == 'normal':
                cellList[i].soma(0.5).gammatoset_bcellcha = gkatp[0] + np.random.normal(0.0,1.0)*np.sqrt(gkatp[1])
            else:
                cellList[i].soma(0.5).gammatoset_bcellcha = gkatp[0]
        def defineBetaHub(cellList,i,hubgkatp=10.0,applytime=5e3):
            # define beta hub cell
            cellList.append(h.betacell())
            setOrigin(cellList[i])
            cellList[i].soma(0.5).gammaapplytime_bcellcha = applytime
            cellList[i].soma(0.5).gammatoset_bcellcha = hubgkatp
    elif model==4:
        ## Created by Chon Lei
        ## The current version is based on the model from
        ## Cha et al. 2011 The Journal of General Physiology and Hraha et al. 2014 PLOS Computational Biology
        ## Last updated: 21/08/2017
        pathToModel = "../models/betacell_cha2011_hraha2014/"
        #TODO
        modelSetup.setupHetDict(varp=hetVar)
        #loadHetero = modelSetup.setHeteroCha2011
        setHetero = modelSetup.setHeteroCha2011
        setOrigin = modelSetup.setOriginCha2011
        HetDict = modelSetup.HetDictCha2011
        def defineBeta(cellList,i,gamma=(0.5,0.0),useDistribution=None,applytime=5e3):
            # define beta cell
            cellList.append(h.betacell())
            cellList[i].soma(0.5).gammaapplytime_bcellcha = applytime
            if hetVar>0:
                setHetero(cellList[i],HetMatrix,i)
            else:
                setOrigin(cellList[i])
            if useDistribution == None:
                cellList[i].soma(0.5).gammamut_bcellcha = gamma[0]
            elif useDistribution == 'sq':
                cellList[i].soma(0.5).gammamut_bcellcha = np.random.uniform(gamma[0],gamma[1])
            elif useDistribution == 'normal':
                cellList[i].soma(0.5).gammamut_bcellcha = gamma[0] + np.random.normal(0.0,1.0)*np.sqrt(gamma[1])
            else:
                cellList[i].soma(0.5).gammamut_bcellcha = gamma[0]
        def defineBetaHub(cellList,i,hubgamma=1.0,applytime=5e3):
            # define beta hub cell
            cellList.append(h.betacell())
            setOrigin(cellList[i])
            cellList[i].soma(0.5).gammaapplytime_bcellcha = applytime
            cellList[i].soma(0.5).gammamut_bcellcha = hubgamma
    elif model==5:
        ## Created by Chon Lei
        ## The current version is based on the model from
        ## Cha et al. 2011 The Journal of General Physiology and Hraha et al. 2014 PLOS Computational Biology
        ## Last updated: 21/08/2017
        pathToModel = "../models/betacell_cha2011_vMetabolic2/"
        #TODO
        modelSetup.setupHetDict(varp=hetVar)
        #loadHetero = modelSetup.setHeteroCha2011
        setHetero = modelSetup.setHeteroCha2011all
        setOrigin = modelSetup.setOriginCha2011all
        HetDict = modelSetup.HetDictCha2011all
        def defineBeta(cellList,i,glu=(0.5,0.0),useDistribution=None,applytime=5e3):
            # define beta cell
            cellList.append(h.betacell())
            cellList[i].soma(0.5).gammaapplytime_bcellcha = applytime
            if hetVar>0:
                setHetero(cellList[i],HetMatrix,i)
            else:
                setOrigin(cellList[i])
            if useDistribution == None:
                cellList[i].soma(0.5).gammatoset_bcellcha = glu[0]
            elif useDistribution == 'sq':
                cellList[i].soma(0.5).gammatoset_bcellcha = np.random.uniform(glu[0],glu[1])
            elif useDistribution == 'normal':
                cellList[i].soma(0.5).gammatoset_bcellcha = glu[0] + np.random.normal(0.0,1.0)*np.sqrt(glu[1])
            else:
                cellList[i].soma(0.5).gammatoset_bcellcha = glu[0]
            #cellList[i].soma(0.5).gKATP_bcellcha = 2.31*1.5 + np.random.normal(0.0,1.0)*2.31*hetVar
            #cellList[i].soma(0.5).PCaL_bcellcha = 48.9*0.7 + np.random.normal(0.0,1.0)*48.9*hetVar
        def defineBetaHub(cellList,i,hubglu=1.0,applytime=5e3):
            # define beta hub cell
            cellList.append(h.betacell())
            setOrigin(cellList[i])
            cellList[i].soma(0.5).gammaapplytime_bcellcha = applytime
            cellList[i].soma(0.5).gammatoset_bcellcha = hubglu
            #cellList[i].soma(0.5).gKATP_bcellcha = 2.31*0.5 + np.random.normal(0.0,1.0)*2.31*hetVar
            #cellList[i].soma(0.5).PCaL_bcellcha = 48.9*1.5 + np.random.normal(0.0,1.0)*48.9*hetVar


    # Delta cell model
    def defineDelta(cellList,i):
        # define delta cell
        cellList.append(h.deltacell())
        #setOrigin(cellList[i])  #TODO: check if we need this


    if gjmodel==1:
        pathToGJModel = "../models/gapjunction_pedersen2015/"
    elif gjmodel==2:
        pathToGJModel = "../models/gapjunction_hermann2010/"

    if species==0:
        pathToMorphology = "../morphologies/mouse/Mouse 40-%d.txt"%morphology
    elif species==1:
        pathToMorphology = "../morphologies/human/H4-1-%d.txt"%morphology
    elif species==2:
        pathToMorphology = "../morphologies/cubic_lattice/cubic%d.txt"%morphology

    random.seed(pyseed)
    np.random.seed(pyseed)
    modelSetup.SetRandomSeed(pyseed)


    ######
    ## Import system setup files (.hoc files and system matrix)
    ######
    CoorDataTmp = np.loadtxt(pathToMorphology)
    # process CoorData to be acceptable format in modelSetup.genCoupleMatrix()
    if mode==3: #species==1 and morphology>3:
        CoorData = CoorDataTmp[np.abs(CoorDataTmp[:,0]-11)<0.1][:,1:4]
        # Load delta cells' coordinates
        DeltaCoorData = CoorDataTmp[np.abs(CoorDataTmp[:,0]-6)<0.1][:,1:4]
        # Control the number of delta cells
        temp = range(len(DeltaCoorData))
        random.shuffle(temp)
        DeltaCoorData = DeltaCoorData[temp[0:50],:]  # say 50 of delta cells
        # Begin: Work out the coupling matrix for delta cells and beta cells
        # (row for delta cell index and column for beta cell index)
        DeltaCoupledMatrix = np.zeros([np.shape(DeltaCoorData)[0], np.shape(CoorData)[0]])
        for i in xrange(np.shape(DeltaCoorData)[0]):
            for j in xrange(np.shape(CoorData)[0]):
                diff = DeltaCoorData[i,:] - CoorData[j,:]
                DeltaCoupledMatrix[i,j] = np.inner(diff,diff)
        HubListByDeltaCell = np.argmin(DeltaCoupledMatrix, axis=1)
        DeltaCoupledMatrix = np.zeros([np.shape(DeltaCoorData)[0], np.shape(CoorData)[0]])
        for i in xrange(np.shape(DeltaCoorData)[0]):
            DeltaCoupledMatrix[i, HubListByDeltaCell[i]] = 1
        # End
        # Return: HubListByDeltaCell, DeltaCoupledMatrix
        ndeltacells = DeltaCoupledMatrix.shape[0]
    else:
        CoorData = CoorDataTmp[CoorDataTmp[:,0]==11][:,1:4] # assume all use `11' as beta cell
    CoupledMatrix = modelSetup.genCoupleMatrix(CoorData,dthres,isletsize,True)
    tempCoupledMatrix = CoupledMatrix + CoupledMatrix.T
    ncells = CoupledMatrix.shape[0]
    nSpatialLinks = np.sum(CoupledMatrix+CoupledMatrix.T,1)
    if ncells != CoupledMatrix.shape[1]:
        raise Exception("CoupledMatrix invalid dimensions.")
    Total = (ncells*ncells)/2 - ncells # maximum number of gapjunctions that could be formed

    if isImitateExp==1:
        imagedCells = modelSetup.getImagedCellIdx(CoorData,topDir=2,imageDepth=10,Ncells=100,method=0)
    else:
        imagedCells = []
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
    try:
        # add delta cell
        h('nrn_load_dll("%sx86_64/.libs/libnrnmech.so")'%'./deltacell/')
        h.load_file (path.join('./deltacell', 'deltacell.hoc'))
    except Exception:
        raise Exception("Please make sure files has been compiled using \n$ nrnivmodl\n")

    ######
    ## System set-up
    ######
    print("*************************")
    print("Starting system set-up...")

    # Define beta hubs cells
    numHubs = int(pHubs*ncells) if pHubs<=1 else pHubs
    if methodToPickHubs!=3:
        temp = range(ncells)
        random.shuffle(temp)
        hubsList = temp[0:numHubs]
        if hubsList_temp!=None and hubsList_temp!=[]:
            hubsList = hubsList_temp
    else: # Pick hubs base on delta cell
        hubsList = list(HubListByDeltaCell)
        if len(HubListByDeltaCell) > numHubs:
            print('Warning: the number of beta hub cells defined by delta ' +
                    '(%d) '%len(HubListByDeltaCell) +
                    'cells is more than required number of hubs (%d)'%numHubs)
        else:
            # Add random beta cell to hub until it reach the required number of hubs
            temp = range(ncells)
            for hub_by_delta in HubListByDeltaCell:
                while temp.count(hub_by_delta) > 0:
                        temp.remove(hub_by_delta)
            random.shuffle(temp)
            hubsList += temp[0:(numHubs-len(HubListByDeltaCell))]
    if False:
        truehub = hubsList[0]
        #truehub2 = hubsList[1]
        hubsList = [truehub]
        hubsList += list(np.arange(ncells)[tempCoupledMatrix[:,truehub]>0])
        #hubsList += [truehub2]
        #hubsList += list(np.arange(ncells)[tempCoupledMatrix[:,truehub2]>0])
    if False:
        # pick #link and dist
        tmpList = []
        print hubsList
        for hubtmp in hubsList:
            if np.sum(tempCoupledMatrix[:,hubtmp])==6:
                tmpList.append(hubtmp)
        hubsList = tmpList[:]
        print hubsList
        for iii in range(len(tmpList)-1):
            for jjj in range(iii+1,len(tmpList)):
                #print np.sqrt(np.sum((CoorData[tmpList[iii]]-CoorData[tmpList[jjj]])**2))
                if np.sqrt(np.sum((CoorData[tmpList[iii]]-CoorData[tmpList[jjj]])**2))>(2.1*dthres) and (iii in hubsList):
                    try:
                        print jjj
                        hubsList.remove(tmpList[jjj])
                    except:
                        pass
        hubsList = hubsList[0:numHubs]
    imagedHubs = list(set(hubsList).intersection(imagedCells))
    numImHubs = int(pHubs*len(imagedCells)) if pHubs<=1 else len(imagedCells)/ncells*pHubs
    if isImitateExp and len(imagedHubs) < max(numImHubs,1) and True:
        nMorehubs = max(numImHubs,1) - len(imagedHubs)
        tempCellsToPick = [x for x in imagedCells if x not in imagedHubs]
        random.shuffle(tempCellsToPick)
        hubsList += tempCellsToPick[0:nMorehubs]
        imagedHubs += tempCellsToPick[0:nMorehubs]
    # Use a previously generated indices
    # Uncomment this to set specific hubs
    #hubsList = [self_def_list]
    #imagedHubs = hubsList[:]
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
    if False:
        aroundHub = []
        for hubhub in hubsList:
            aroundHub += list(np.arange(ncells)[tempCoupledMatrix[:,hubhub]>0])
        temp = [i for i in temp if i not in aroundHub]
    random.shuffle(temp)
    nonHubsToPickList = temp[0:numNonHubsToPick]
    imagedNonHubs = []
    if isImitateExp:
        temp = [i for i in imagedCells if i not in imagedHubs]
        random.shuffle(temp)
        imagedNonHubs = temp[0:len(imagedHubs)]
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
    def silenceCellV(iclampList,cell,delay=250e3,dur=250e3,amp=-100.0):
        # vclamp cell to imitate cell silencing in experiments
        iclampList.append(h.SEClamp (0.5, sec = cell.soma) )
        iclampList[-1].dur1 = silenceStart
        iclampList[-1].rs = 1e30
        iclampList[-1].dur2 = dur
        iclampList[-1].amp2 = amp
    
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
    toPick = random.randint(0,len(hubsList))
    for i in range(ncells):
        if i not in hubsList:
            defineBeta(cell,i,**(model_kwargs['beta']))
            if isImitateExp:
                #if i in list(np.arange(ncells)[tempCoupledMatrix[:,imagedHubs[whichHub]]>0]):
                #if (i == imagedNonHubs[whichHub]) and (mode==2):
                #if False and (i in imagedCells[:50]):
                if mode==2 and i in nonHubsToPickList[:]:
                    print "silencing cell ",i
                    with open(outlog, 'a') as f:
                        f.write('#silencedCell = %d\n'%i)
                        f.write('#cell%d_nSpatialLinks = %d\n'%(i,nSpatialLinks[i]))
                    silenceCell(iclamp_hubs,cell[i],silenceStart,silenceDur,silenceAmp)
                    raise Exception('This simulation should not go in here...')
            else:
                if (i in nonHubsToPickList[:]) and (mode==2):
                    print "silencing cell ",i
                    with open(outlog, 'a') as f:
                        f.write('#silencedCell = %d\n'%i)
                        f.write('#cell%d_nSpatialLinks = %d\n'%(i,nSpatialLinks[i]))
                    silenceCell(iclamp_hubs,cell[i],silenceStart,silenceDur,silenceAmp)
                    raise Exception('This simulation should not go in here...')
        else:
            defineBetaHub(cell,i,**(model_kwargs['betahub']))
            if isImitateExp:
                if mode==1 and (i in hubsList[:4]):#i==imagedHubs[whichHub]:
                    #or i in list(np.arange(ncells)[tempCoupledMatrix[:,imagedHubs[whichHub]]>0]):
                    print "silencing cell ",i
                    with open(outlog, 'a') as f:
                        f.write('#silencedCell = %d\n'%i)
                        f.write('#cell%d_nSpatialLinks = %d\n'%(i,nSpatialLinks[i]))
                    silenceCell(iclamp_hubs,cell[i],silenceStart,silenceDur,silenceAmp)
                    raise Exception('This simulation should not go in here...')
            else:
                if mode==1 and i==hubsList[toPick]:
                    print "silencing cell ",i
                    with open(outlog, 'a') as f:
                        f.write('#silencedCell = %d\n'%i)
                        f.write('#cell%d_nSpatialLinks = %d\n'%(i,nSpatialLinks[i]))
                    silenceCell(iclamp_hubs,cell[i],silenceStart,silenceDur,silenceAmp)
                    raise Exception('This simulation should not go in here...')


    # Define delta cell
    if mode==3:
        print "Defining delta cells..."
        assert iclamp_hubs == []  # Make sure we are doing the right thing...
        delta_cell = []
        syn = []
        for i in range(ndeltacells):
            defineDelta(delta_cell,i)
            # Begin: Set connection between the neighbour hub
            syn.append(h.sst(0.5, sec=cell[HubListByDeltaCell[i]].soma))
            h.setpointer(delta_cell[i].soma(0.5)._ref_T_rel, 'cp', syn[i])
            # End
            print "Injecting current (IClamp) to delta cell ",i
            with open(outlog, 'a') as f:
                f.write('#IClampingDeltaCell = %d\n'%i)
            silenceCell(iclamp_hubs,delta_cell[i],silenceStart,silenceDur,silenceAmp)
    

    # Declare heterogeneity GJ matrix
    HetGjMatrix = np.zeros(CoupledMatrix.shape)
    # Use a previously generated heterogeneity GJ matrix
    #HetMatrix = np.loadtxt('')

    
    #TODO can put this in modelSetup.py too
    print "Defining gap junction connections..."
    gap = []
    for i in range(ncells):
        for j in range(ncells):
            a = np.random.rand()
            if CoupledMatrix[i,j] > 0 and ((i in hubsList) or (j in hubsList)) and a<p_connect:
                if pggaphubstd > 0:
                    HetGjMatrix[i,j] = max(ggaphub * (1.0 + np.random.normal(0.0,1.0)*pggaphubstd), 0)
                else:
                    HetGjMatrix[i,j] = ggaphub
                if False:
                    if i in hubsList:
                        HetGjMatrix[i,j] = ggaphub*7.0/np.sum(tempCoupledMatrix[:,i]>0) * (1.0 + np.random.normal(0.0,1.0)*pggaphubstd)
                    else:
                        HetGjMatrix[i,j] = ggaphub*7.0/np.sum(tempCoupledMatrix[:,j]>0) * (1.0 + np.random.normal(0.0,1.0)*pggaphubstd)
                gap.append(h.gapjunction(cell[i], cell[j], 0.5, 0.5, HetGjMatrix[i,j]*CoupledMatrix[i,j],gjtau))
            elif CoupledMatrix[i,j] > 0 and a<p_connect: #and ((i not in nonHubsToPickList) or (j not in nonHubsToPickList)):
                if pggapstd > 0:
                    HetGjMatrix[i,j] = max(ggap * (1.0 + np.random.normal(0.0,1.0)*pggapstd), 0)
                else:
                    HetGjMatrix[i,j] = ggap
                gap.append(h.gapjunction(cell[i], cell[j], 0.5, 0.5, HetGjMatrix[i,j]*CoupledMatrix[i,j],gjtau))


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
    deltavrec = []
    deltacarec = []
    for i in range(ncells):
        vrec.append(h.Vector())
        carec.append(h.Vector())
        vrec[i].record(cell[i].soma(0.5)._ref_v)
        carec[i].record(cell[i].soma(0.5)._ref_cai)
    for i in range(ndeltacells):
        deltavrec.append(h.Vector())
        deltacarec.append(h.Vector())
        deltavrec[i].record(delta_cell[i].soma(0.5)._ref_v)
        deltacarec[i].record(delta_cell[i].soma(0.5)._ref_cai)


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
        if temptstop >= silenceStart and False:
           for iclamp in iclamp_hubs:
                iclamp.rs = 0.001
        #if temptstop > 150e3:
        #    for iclamp in iclamp_hubs:
        #        iclamp.amp = -10.0
        #if temptstop > 200e3:
        #    for iclamp in iclamp_hubs:
        #        iclamp.amp = -20.0
        #if temptstop > 250e3:
        #    for iclamp in iclamp_hubs:
        #        iclamp.amp = -30.0
        #if temptstop > 300e3:
        #    for iclamp in iclamp_hubs:
        #        iclamp.amp = -40.0
        temptstop += tbatch  # tstop for current batch
        h.frecord_init()  # reuse all recording vectors
        h.continuerun(temptstop)
        # exporting Ca time series
        tosave = modelSetup.convertSimOutput(carec,downSampling,reuse=True)
        modelSetup.savedat(outCa,tosave,'Ca',outlog,idx=i)
        # exporting Vm time series
        tosave = modelSetup.convertSimOutput(vrec,downSampling,reuse=True)
        modelSetup.savedat(outVm,tosave,'Vm',outlog,idx=i)
        # exporting Vm time series
        tosave = modelSetup.convertSimOutput(deltacarec,downSampling,reuse=True)
        modelSetup.savedat(outDeltaCa,tosave,'DeltaCa',outlog,idx=i)
        # exporting Vm time series
        tosave = modelSetup.convertSimOutput(deltavrec,downSampling,reuse=True)
        modelSetup.savedat(outDeltaVm,tosave,'DeltaVm',outlog,idx=i)
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
        # exporting Vm time series
        tosave = modelSetup.convertSimOutput(deltacarec,downSampling,reuse=True)
        modelSetup.savedat(outDeltaCa,tosave,'DeltaCa',outlog,idx=i+1)
        # exporting Vm time series
        tosave = modelSetup.convertSimOutput(deltavrec,downSampling,reuse=True)
        modelSetup.savedat(outDeltaVm,tosave,'DeltaVm',outlog,idx=i+1)
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


if __name__== '__main__':
    try:
        print("Loading external model parameters from %s"%sys.argv[1])
        modelParam = pickle.load(open(sys.argv[1],'rb'))
    except Exception:
        print("Cannot read external model parameters...\nUsing default parameters...\n\n")
    #modelParam['pHubs'] = 0.05
    #modelParam['whichHub'] = 3
    #modelParam['ggaphub'] = 0.4
    #modelParam['silenceAmp'] = -30.0
    # just try to run main(), but not expecting to be used in this way
    main(modelParam)

## 
