import numpy as np
import itertools
import sys

try:
    idx = int(sys.argv[1])
    N = int(sys.argv[2])
    dist = float(sys.argv[3])
except Exception:
    idx = 0
    N = 5
    dist = 17.0

toIter = [xrange(N) for i in xrange(3)]
lattice = np.array(list(itertools.product(*toIter)))*dist

thefile = open('cubic%d.txt'%idx, 'w')
for i in xrange(len(lattice)):
    thefile.write("11")
    for a in lattice[i]:
        thefile.write("\t%f"%a)
    thefile.write("\n")

# end
