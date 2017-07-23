from neuron import h
import numpy as np
import matplotlib.pylab as plt

## 
## Simulation of the electrophysiology of the beta-cell.
## Using package NEURON yale.
## 
## Created by Chon Lei
## The currently version is based on the model from 
## Morten Gram Pedersen 2010 Biophysical Journal Vol. 99
## Last updated: 10/02/2017
## 


## Import system setup files (.hoc files and system matrix)
CoupledMatrix = np.loadtxt('CouplingMatrix.dat',delimiter=' ')
ncells = CoupledMatrix.shape[0]
if ncells != CoupledMatrix.shape[1]:
    raise Exception("CoupledMatrix invalid dimensions.")

h.load_file ("betacell.hoc")
h.load_file("gapjunction.hoc")

Total = (ncells*ncells)/2 - ncells # maximum number of gapjunctions that could be formed


## System set-up
print("*************************")
print("Starting system set-up...")
cell = []
for i in range(ncells):
    cell.append(h.betacell())

gap = []
for i in range(ncells):
    for j in range(ncells):
        if CoupledMatrix[i,j] > 0:
            gap.append(h.gapjunction(cell[i], cell[j], 0.5, 0.5, 0.00017e-1*CoupledMatrix[i,j]))  # 170pS (Zhang et al. 2008) Why the factor 1e-4


## External stimulation
"""
stimulus = h.IClamp (0.5, sec = cell[0].soma)
stimulus.delay = 100
stimulus.dur = 300
stimulus.amp = 0.5
"""


## System recorder initialisation
time = h.Vector()
time.record(h._ref_t)
vrec = []
for i in range(ncells):
    vrec.append(h.Vector())
    vrec[i].record(cell[i].soma(0.5)._ref_v)


## Main simulation
h.load_file("stdrun.hoc")
h.init()
# h.v_init = purkinjecelly.undershootpotential
h.tstop = 1e3
h.dt = 0.1
h.steps_per_ms = 1./h.dt
print("Starting main simulation...")
h.run()
print("Simulation completed! :)")
print("*************************")


## Exporting


## Visualisation
time = np.array(time)
membranepotential = np.array(vrec[0])
plt.plot(time, membranepotential)

membranepotential2 = np.array(vrec[102])
plt.plot(time, membranepotential2)

plt.show()
