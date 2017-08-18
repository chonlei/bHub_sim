from neuron import h
import numpy as np
import matplotlib.pylab as plt

## 
## Simulation of the electrophysiology of the beta-cell.
## Using package NEURON yale.
## 
## Created by Chon Lei
## The currently version is based on the model from 
## Cha et al. 2011
## Last updated: 25/03/2017
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

vm1 = h.Vector()
vm1.record (betacell.soma(0.5)._ref_v)
vm2 = h.Vector()
vm2.record (betacell2.soma(0.5)._ref_v)
"""
atp = h.Vector()
atp.record (betacell.soma(0.5)._ref_atpi)
adp = h.Vector()
adp.record (betacell.soma(0.5)._ref_mgadpi)

"""
ca = h.Vector()
ca.record (betacell.soma(0.5)._ref_cai)

na = h.Vector()
na.record (betacell.soma(0.5)._ref_nai)
k = h.Vector()
k.record (betacell.soma(0.5)._ref_ki)

## Main simulation
h.load_file("stdrun.hoc")
h.init()
h.v_init = -48.9045  #-69.8663703359279 or #-48.9045
h.tstop = 25e4
h.dt = 0.1
h.steps_per_ms = 1./h.dt
print("Starting main simulation...")
h.run()
print("Simulation completed! :)")
print("*************************")


## Exporting


## Visualisation
time = np.array(time)

plt.figure(1)
vm1 = np.array(vm1)
plt.plot(time, vm1, 'r-',label='cell1')
vm2 = np.array(vm2)
plt.plot(time, vm2, 'b-',label='cell2')
plt.legend()
plt.title("V", fontsize=20)
plt.xlabel("time [ms]", fontsize=20)
plt.ylabel("V [mV]", fontsize=20)

"""
plt.figure(2)
atp = np.array(atp)
plt.plot(time, atp, 'b-',label='atp')
plt.legend()
plt.title("ATP", fontsize=20)
plt.xlabel("time [ms]", fontsize=20)
plt.ylabel("conc [mM]", fontsize=20)


plt.figure(3)
adp = np.array(adp)
plt.plot(time, adp, 'r-',label='adp')
plt.legend()
plt.title("MgADP", fontsize=20)
plt.xlabel("time [ms]", fontsize=20)
plt.ylabel("conc [mM]", fontsize=20)

"""
plt.figure(4)
ca = np.array(ca)
plt.plot(time, ca, 'r-',label='Ca')
plt.legend()
plt.title("Ca", fontsize=20)
plt.xlabel("time [ms]", fontsize=20)
plt.ylabel("conc [mM]", fontsize=20)


plt.figure(5)
na = np.array(na)
plt.plot(time, na, 'b-',label='Na')
#k = np.array(k)
#plt.plot(time, k, 'k-',label='K')
plt.legend()
plt.title("Na", fontsize=20)
plt.xlabel("time [ms]", fontsize=20)
plt.ylabel("conc [mM]", fontsize=20)


plt.show()