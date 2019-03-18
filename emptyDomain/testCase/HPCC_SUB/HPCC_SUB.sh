#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N LES
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe mpi 36
#$ -P quanah

JOBDIR='/lustre/work/wan39502/emptyDomain'

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Mesh generation
# cd ${JOBDIR}/testCase
# runApplication blockMesh
# runApplication checkMesh

# Mesh info extraction
# wmake ${JOBDIR}/meshInfo
# ${JOBDIR}/bin/meshInfo

# Inflow turbulence generation
# cd ${JOBDIR}/turbGen
# module load matlab/R2018b
# matlab -nodisplay -nosplash < turbGen.m | tee log.turbGen

# Inflow mass flux correction
# cd ${JOBDIR}/massFluxCorrection
# python massFluxCorrection.py | tee log.massFluxCorrection

cd ${JOBDIR}/testCase

[ ! -d 0 ] && cp -r 0.orig 0

# Decompose the mesh
# runApplication decomposePar

# pisoFoam solver
runParallel $(getApplication)

# Reconstructs fields of a case that is decomposed for parallel execution
runApplication reconstructPar -latestTime
