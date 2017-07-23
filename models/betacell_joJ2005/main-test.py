from neuron import h
import numpy as np
import matplotlib.pylab as plt


arr = np.loadtxt('CouplingMatrix.dat',delimiter=' ')


h.load_file ("morphology.hoc")
h.load_file("gapjunction.hoc")

betacell = h.CellTemp()
betacell2 = h.CellTemp()
#g = h.gapjunction(betacell,betacell2,0.5,0.5,0.00017e-1)


stimulus = h.IClamp (0.5, sec = betacell.soma)
stimulus.delay = 100
stimulus.dur = 300
stimulus.amp = 0.1   



time = h.Vector()
membranepotential = h.Vector()
time.record(h._ref_t)
membranepotential.record (betacell.soma(0.5)._ref_v)
membranepotential2 = h.Vector()
membranepotential2.record (betacell2.soma(0.5)._ref_v)



h.load_file("stdrun.hoc")
h.init()
# h.v_init = purkinjecelly.undershootpotential
h.tstop = 500
h.run()



time = np.array(time)
membranepotential = np.array(membranepotential)
plt.plot(time, membranepotential)

membranepotential2 = np.array(membranepotential2)
plt.plot(time, membranepotential2)

plt.show()
