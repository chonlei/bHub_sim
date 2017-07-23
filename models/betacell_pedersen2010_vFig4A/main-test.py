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


## Import system setup files (.hoc files)
h.load_file ("betacell.hoc")
h.load_file("gapjunction.hoc")


## System set-up
print("*************************")
print("Testing system: two coupled beta cells.")
print("Starting system set-up...")
betacell = h.betacell()
betacell2 = h.betacell()
#betacell2.soma.gbar_cal = 0.14e-4
g = h.gapjunction(betacell,betacell2,0.5,0.5,0.00017e-1)


## External stimulation
'''
stimulus = h.IClamp (0.5, sec = betacell.soma)
stimulus.delay = 100
stimulus.dur = 300
stimulus.amp = 0.1   
'''


## System recorder initialisation
time = h.Vector()
time.record(h._ref_t)
membranepotential = h.Vector()
membranepotential.record (betacell.soma(0.5)._ref_v)
membranepotential2 = h.Vector()
membranepotential2.record (betacell2.soma(0.5)._ref_v)


## Main simulation
h.load_file("stdrun.hoc")
h.init()
# h.v_init = purkinjecelly.undershootpotential
h.tstop = 1e4
h.dt = 0.1
h.steps_per_ms = 1./h.dt
print("Starting main simulation...")
h.run()
print("Simulation completed! :)")
print("*************************")


## Exporting


## Visualisation
time = np.array(time)
membranepotential = np.array(membranepotential)
plt.plot(time, membranepotential)

membranepotential2 = np.array(membranepotential2)
plt.plot(time, membranepotential2)
plt.xlabel("t [ms]")
plt.ylabel("V [mV]")

plt.show()
