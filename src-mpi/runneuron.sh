#!/bin/bash

#SBATCH --job-name=neuron-mpi
#SBATCH --time=10:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --partition=compute

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=linfordbriant@gmail.com 


module load neuron

# create mpirun argument
HOST_LIST=$(scontrol show hostnames $SLURM_JOB_NODELIST | tr "\n" "," | sed 's/.$//')
MPI_HOSTS=" -np "$SLURM_NTASKS" -hosts "$HOST_LIST""

# change directory to source code
cd ${HOME}/bHub_sim/src-mpi

# create output directory
$OUTPUTINDEX=1
outputdir=$(python outputSetup.py ${DATA}/output/ ${OUTPUTINDEX} >&1)
# run main simulations
$SEED=1
$MODE=1
$ISLET=1
$GJMODEL=2
$GJHUB=0.05
$GJNONHUB=0.01
$SPECIES=0
$ISACT=0
$NHUB=3
mpirun $MPI_HOSTS python main.py ${outputdir} $SEED $MODE $ISLET $GJMODEL $GJHUB $GJNONHUB $SPECIES $ISACT $NHUB
