#!/bin/bash
# --------------------------------------------------------------------------- #
#   ==   ===  ==                                                              #
#  ||   ||=  || ))  support s. r. o. 2017, www.cfdsupport.com                 #
#   ==        ==                                                              #
# --------------------------------------------------------------------------- #


if [[ "$1" -eq "" ]] ; then
    numProcs=6
else
    numProcs=$1
fi

echo numProcs
echo $numProcs

sed -i "s/\(.*numberOfSubdomains[ \t]*\)[0-9].*;/numberOfSubdomains $numProcs;/g" mesh/nonRotatingPart/system/decomposeParDict
sed -i "s/\(.*numberOfSubdomains[ \t]*\)[0-9].*;/numberOfSubdomains $numProcs;/g" mesh/rotatingPart/system/decomposeParDict
sed -i "s/\(.*n[ \t]*\)(.*;/\1(1 1 $numProcs);/g" mesh/nonRotatingPart/system/decomposeParDict

echo "Creating mesh of rotating part..."
   cd  mesh/rotatingPart
   blockMesh > ../../log.rotatingPart-blockMesh  2>&1
   surfaceFeatureExtract > ../../log.rotatingPart-surfaceFeatureExtract 2>&1
   decomposePar > ../../log.rotatingPart-decomposePar  2>&1
   mpiexec -np $numProcs transformPoints -parallel -rotate '( (1 0 0) (0.993572 0 0.113203) )'> ../../log.rotatingPart-transformPoints  2>&1
   mpiexec -np $numProcs snappyHexMesh -overwrite -parallel > ../../log.rotatingPart-snappyHexMesh  2>&1
   reconstructParMesh -constant > ../../log.rotatingPart-reconstructParMesh  2>&1
   topoSet > ../../log.rotatingPart-topoSet 2>&1 # creating 'rotor' zone
   checkMesh  -constant > ../../log.rotatingPart-checkMesh  2>&1
   cd ..

echo "Creating mesh of non rotating part..."
   cd  nonRotatingPart
   blockMesh > ../../log.nonRotatingPart-blockMesh  2>&1
   surfaceFeatureExtract > ../../log.nonRotatingPart-surfaceFeatureExtract 2>&1
   decomposePar > ../../log.nonRotatingPart-decomposePar  2>&1
   mpiexec -np $numProcs snappyHexMesh -overwrite -parallel > ../../log.nonRotatingPart-snappyHexMesh  2>&1
   reconstructParMesh -constant > ../../log.nonRotatingPart-reconstructParMesh  2>&1
   topoSet > ../../log.nonRotatingPart-topoSet 2>&1                                                              # creating 'Co1_zone' zone
   checkMesh  -constant > ../../log.nonRotatingPart-checkMesh  2>&1
   cd ..

echo "Merging meshes..."
   mkdir final
   cp -rv nonRotatingPart/constant final/ > ../log.finalMesh-cp 2>&1
   cp -rv nonRotatingPart/system final/ >> ../log.finalMesh-cp 2>&1
   mergeMeshes final rotatingPart -overwrite > ../log.finalMesh-mergeMeshes 2>&1
   cd final
   createPatch -overwrite > ../../log.finalMesh-createPatch 2>&1
   changeDictionary -dict "../../system/changeDictionaryDict" > ../../log.caseMesh-changeDictionary 2>&1
   transformPoints -scale '(0.001 0.001 0.001)' > ../../log.finalMesh-transformPoints 2>&1
   cd ../..

echo "Copying mesh into case directory..."
   cp -rv mesh/final/constant/polyMesh constant > log.finalMesh-cpMeshToCase 2>&1

