#!/bin/bash
# --------------------------------------------------------------------------- #
#   ==   ===  ==                                                              #
#  ||   ||=  || ))  support s. r. o. 2017, www.cfdsupport.com                 #
#   ==        ==                                                              #
# --------------------------------------------------------------------------- #
cd ${0%/*} || exit 1    # run from this directory

echo "--------"
echo "Cleaning tutorials ..."
find . \( -name 'processor[0-9]*' \) -exec rm -rfv {} \;
rm -rfv log.*
rm -rfv mesh/final
rm -rfv mesh/nonRotatingPart/constant/polyMesh
rm -rfv mesh/nonRotatingPart/constant/extendedFeatureEdgeMesh
rm -rfv mesh/nonRotatingPart/constant/triSurface/*.eMesh
rm -rfv mesh/rotatingPart/constant/polyMesh
rm -rfv mesh/rotatingPart/constant/extendedFeatureEdgeMesh
rm -rfv mesh/rotatingPart/constant/triSurface/*.eMesh
rm -rfv postProcessing
rm -rfv 0
rm -rfv constant/polyMesh
echo "--------"

# ----------------------------------------------------------------- end-of-file
