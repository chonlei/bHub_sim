#!/usr/bin/env python

import numpy as np
import os
import glob
from mpi4py import MPI
from neuron import h
import simulator

'''
pc = h.ParallelContext()
id = int(pc.id())
nhost = int(pc.nhost())
print "I am ", id, " of ", nhost
'''

def main(comm):
    # Get rank and size
    size = comm.Get_size()
    rank = comm.Get_rank()
    # check they are doing right thing
    #print("RANK %d of SIZE %d is doing the right job..."%(rank,size))
    simulator.main()


if __name__ == "__main__":
    main(MPI.COMM_WORLD)
