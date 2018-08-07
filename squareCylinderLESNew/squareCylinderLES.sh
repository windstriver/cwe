#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N LES_NEW3
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe mpi 144
#$ -P quanah

cd testCase

#source $HOME/OpenFOAM/OpenFOAM-5.0/etc/bashrc

#blockMesh

#checkMesh

#../bin/writePatchCentres

decomposePar

mpirun -np $NSLOTS pimpleFoam -parallel

#reconstructPar

#foamToEnsight -fields '("p.*" "U.*")'

