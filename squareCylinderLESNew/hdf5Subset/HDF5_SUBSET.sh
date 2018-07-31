#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N HDF5_SUBSET
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe mpi 36
#$ -P quanah

module load python3/3.6.4

HDF5_USE_FILE_LOCKING=FALSE

python hdf5Subset.py
