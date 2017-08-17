#!/bin/bash

#SBATCH --job-name=neuron-mpi
#SBATCH --time=100:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=16
#SBATCH --partition=compute

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=chon.lei@cs.ox.ac.uk 


module load neuron

# create mpirun argument
HOST_LIST=$(scontrol show hostnames $SLURM_JOB_NODELIST | tr "\n" "," | sed 's/.$//')
MPI_HOSTS=" -np "$SLURM_NTASKS" -hosts "$HOST_LIST""

# change directory to source code
cd ${HOME}/bHub_sim/src-mpi

# create output directory
outputdir=$(python outputSetup.py ${DATA} >&1)
# run main simulations
mpirun $MPI_HOSTS python main.py ${outputdir}

