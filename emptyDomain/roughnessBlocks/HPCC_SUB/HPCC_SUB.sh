#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N MESH
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe mpi 36
#$ -P quanah

cd $WORK/roughnessBlocksCoarse/

./Allrun
