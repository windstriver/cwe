#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N CDRFG
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe mpi 36
#$ -P quanah

JOBDIR='/lustre/work/wan39502/emptyDomain'

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Mesh generation
#cd ${JOBDIR}/testCase
#runApplication blockMesh
#runApplication checkMesh

# Mesh info extraction
#wmake ${JOBDIR}/writeMesh
##${JOBDIR}/bin/writeMesh

# Inflow turbulence generation
#cd ${JOBDIR}/turbGen
#module load matlab/R2018b
#matlab -nodisplay -nosplash < turbGen.m | tee log.turbGen

# Inflow mass flux correction
#cd ${JOBDIR}/massFluxCorrection
#python massFluxCorrection.py | tee log.massFluxCorrection

cd ${JOBDIR}/testCase

../bin/writeHDF5
