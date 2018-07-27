/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "Create mesh, no clear-out\n" << endl;
    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Return a list of patch names
    // wordList patchNames = mesh.boundaryMesh().names();
    // Info<< "List of patch names: " << patchNames << endl;

    // Find patch index given a patch name
    label patchInlet = mesh.boundaryMesh().findPatchID("inlet");
    label patchBuilding = mesh.boundaryMesh().findPatchID("building");

    const vectorField& faceCentresInlet = mesh.Cf().boundaryField()[patchInlet];
    const vectorField& faceCentresBuilding = mesh.Cf().boundaryField()[patchBuilding];

    // Write patch centres to a file
    fileName outputFile("inletPatchFaceCentres");
    fileName outputFile2("buildingPatchFaceCentres");
    
    // Write inlet face centres to file
    OFstream os(outputFile);
    forAll(faceCentresInlet, faceI)
    {
      os << faceI << ","
    	 << faceCentresInlet[faceI].x() << ","
     	 << faceCentresInlet[faceI].y() << ","
     	 << faceCentresInlet[faceI].z() << endl;
    }

    // Write building face centres for post-processing
    OFstream os2(outputFile2);
    forAll(faceCentresBuilding, faceI)
    {
        os2 << faceCentresBuilding[faceI] << endl;
    }
  
    // Return face area vector as surfaceVectorField
    const vectorField& faceAreaVectors = mesh.Sf().boundaryField()[patchBuilding];
    // Write face area vectors to a file
    fileName outputFile3("buildingPatchFaceAreaVectors");
    OFstream os3(outputFile3);

    forAll(faceAreaVectors, faceI)
    {
        os3 << faceAreaVectors[faceI].x() << "," 
            << faceAreaVectors[faceI].y() << "," 
            << faceAreaVectors[faceI].z() << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
