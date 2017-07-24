from neuron import h
import numpy as np
import matplotlib.pylab as plt

## 
## Simulation of the electrophysiology of the beta-cell.
## Using package NEURON yale.
## 
## Created by Chon Lei
## The currently version is based on the model from 
## M. Meyer-Hermann 2007 Biophysical Journal Vol. 93
## Last updated: 19/02/2017
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
membranepotential = h.Vector()
membranepotential.record (betacell.soma(0.5)._ref_m_skca)
membranepotential2 = h.Vector()
membranepotential2.record (betacell.soma(0.5)._ref_h_pmca)
membranepotential3 = h.Vector()
membranepotential3.record (betacell.soma(0.5)._ref_g_cal)
membranepotential4 = h.Vector()
membranepotential4.record (betacell.soma(0.5)._ref_m_katp)
membranepotential41 = h.Vector()
membranepotential41.record (betacell.soma(0.5)._ref_h_ncx)
membranepotential42 = h.Vector()
membranepotential42.record (betacell.soma(0.5)._ref_g_cat)

membranepotential5 = h.Vector()
membranepotential5.record (betacell.soma(0.5)._ref_cai)

membranepotential6 = h.Vector()
membranepotential6.record (betacell.soma(0.5)._ref_v)
membranepotential7 = h.Vector()
membranepotential7.record (betacell2.soma(0.5)._ref_v)

## Main simulation
h.load_file("stdrun.hoc")
h.init()
h.v_init = -70
h.tstop = 9e4
h.dt = 0.1
h.steps_per_ms = 1./h.dt
print("Starting main simulation...")
h.run()
print("Simulation completed! :)")
print("*************************")


## Exporting
membranepotential6 = np.array(membranepotential6)
np.savetxt("testoutput.txt", membranepotential6)


"""
## Visualisation
time = np.array(time)

plt.figure(1)
membranepotential = np.array(membranepotential)
plt.plot(time, membranepotential, 'g--',label="m(sK,Ca)")

membranepotential2 = np.array(membranepotential2)
plt.plot(time, membranepotential2, 'b-',label="H(PMCA)")

membranepotential3 = np.array(membranepotential3)
plt.plot(time, membranepotential3, 'r-',label="g*(1-H)(Ca,L)")

membranepotential4 = np.array(membranepotential4)
plt.plot(time, (1 - membranepotential4), 'k-.',label="1-m(K,ATP)")

membranepotential41 = np.array(membranepotential41)
plt.plot(time, membranepotential41, 'm--',label="H(NCX)")

membranepotential42 = np.array(membranepotential42)
plt.plot(time, membranepotential42, 'c-',label="g*h(Ca,T)")
plt.legend()
plt.title("Reproduce Hermann 2007 Fig.4A", fontsize=20)
plt.xlabel("time [ms]", fontsize=20)
plt.ylabel("open/action probability", fontsize=20)


plt.figure(2)
membranepotential5 = np.array(membranepotential5)
plt.plot(time, membranepotential5/membranepotential5[0], 'r-')
plt.title("C/C_0", fontsize=20)
plt.xlabel("time [ms]", fontsize=20)
plt.ylabel("C/C_0", fontsize=20)


plt.figure(3)
membranepotential6 = np.array(membranepotential6)
plt.plot(time, membranepotential6, 'b-',label='cell1')
membranepotential7 = np.array(membranepotential7)
plt.plot(time, membranepotential7, 'r-',label='cell2')
plt.legend()
plt.title("V", fontsize=20)
plt.xlabel("time [ms]", fontsize=20)
plt.ylabel("V [mV]", fontsize=20)


plt.show()
"""