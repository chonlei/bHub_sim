from neuron import h
import numpy as np
import matplotlib.pylab as plt

## exactly reproduce betacell from Linford 08/02/2017

CoupledMatrix = np.loadtxt('CouplingMatrix.dat',delimiter=' ')
ncells = CoupledMatrix.shape[0]

h.load_file ("morphology.hoc")
h.load_file("gapjunction.hoc")

Total = (ncells*ncells)/2 - ncells # maximum number of gapjunctions that could be formed



cell = []
for i in range(ncells):
    cell.append(h.CellTemp())



gap = []
for i in range(ncells):
    for j in range(ncells):
        if CoupledMatrix[i,j] > 0:
            gap.append(h.gapjunction(cell[i], cell[j], 0.5, 0.5, 0.00017e-1*CoupledMatrix[i,j]))  # 170pS (Zhang et al. 2008) Why the factor 1e-4



stimulus = h.IClamp (0.5, sec = cell[0].soma)
stimulus.delay = 100
stimulus.dur = 300
stimulus.amp = 0.5



time = h.Vector()
time.record(h._ref_t)
vrec = []
for i in range(ncells):
    vrec.append(h.Vector())
    vrec[i].record(cell[i].soma(0.5)._ref_v)



h.load_file("stdrun.hoc")
h.init()
# h.v_init = purkinjecelly.undershootpotential
h.tstop = 1e3
h.dt = 0.1
h.steps_per_ms = 1./h.dt
h.run()



time = np.array(time)
membranepotential = np.array(vrec[0])
plt.plot(time, membranepotential)

membranepotential2 = np.array(vrec[102])
plt.plot(time, membranepotential2)

plt.show()
