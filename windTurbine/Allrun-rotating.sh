#!/bin/bash
# --------------------------------------------------------------------------- #
#   ==   ===  ==                                                              #
#  ||   ||=  || ))  support s. r. o. 2017, www.cfdsupport.com                 #
#   ==        ==                                                              #
# --------------------------------------------------------------------------- #

# number of CPUs to run on
numProcs=6
# endTime
endTime=1

# check environment
if [[ $(echo $WM_PROJECT_VERSION |  cut -c1-3) != "dev" ]];
then
    echo "Use OpenFOAM dev with this example script, please."
    exit
fi

echo
echo "Cleaning..."
   ./Allclean.sh

sed -i "s/\(.*numberOfSubdomains[ \t]*\)[0-9].*;/numberOfSubdomains $numProcs;/g" system/decomposeParDict
sed -i "s/\(.*endTime[ \t]*\)[0-9].*;/\1$endTime;/g" system/controlDict


   ./makeMesh.sh $numProcs

echo "Running Simulation..."
   cp -r 0.org 0
   decomposePar > log.simulation-decomposePar  2>&1 
   mpiexec -np $numProcs renumberMesh -overwrite -parallel > log.simulation-renumberMesh 2>&1
   mpiexec -np $numProcs pimpleDyMFoam -parallel > log.simulation-pimpleDyMFoam 2>&1
