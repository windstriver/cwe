#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N UI0
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe mpi 36
#$ -P quanah

JOBDIR='/lustre/work/wan39502/bridgeSectionUniformInflow'
CASE='testCase'

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Mesh generation
cd ${JOBDIR}/${CASE}
runApplication blockMesh
runApplication checkMesh
# runApplication snappyHexMesh

# Simulation
cd ${JOBDIR}/${CASE}

[ ! -d 0 ] && cp -r 0.orig 0

# Map fileds for init
# mapFields ../testCase3 -sourceTime "latestTime"

# Decompose the mesh
runApplication decomposePar

# pisoFoam solver
# runApplication $(getApplication)
runParallel $(getApplication)

# Reconstructs fields of a case that is decomposed for parallel execution
runApplication reconstructPar -latestTime

