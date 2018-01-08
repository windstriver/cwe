#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N LES_INLET
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe mpi 36
#$ -P quanah

cd testCase

#source $HOME/OpenFOAM/OpenFOAM-5.0/etc/bashrc

#blockMesh

#../bin/createHDF5

#../bin/writeHDF5

decomposePar

#module load gnu openmpi

mpirun -np $NSLOTS pimpleFoam -parallel

reconstructPar

