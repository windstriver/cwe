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

    // Info<< mesh.C() << endl;
    // Info<< mesh.V() << endl;

    // Return a list of patch names
    // wordList patchNames = mesh.boundaryMesh().names();
    // Info<< "List of patch names: " << patchNames << endl;

    // Find patch index given a patch name
    label patchI = mesh.boundaryMesh().findPatchID("building");
    //Info<< "Patch index for building is: " << patchI << endl;

    // Return face centres as surfaceVectorField
    // surfaceVectorField faceCentres = mesh.Cf();
    // Info<< "surfaceVectorField: " << faceCentres << endl;

    const vectorField& faceCentres = mesh.Cf().boundaryField()[patchI];
    // Info<< "faceCentres on patch fixedWalls" << nl << faceCentres << endl;

    // Write patch centres to a file
    fileName outputFile("buildingPatchFaceCentres");
    // fileName outputFile("inletPatchFaceCentres");
    OFstream os(outputFile);

    forAll(faceCentres, faceI)
    {
        os << faceCentres[faceI] << endl;
    }
    // forAll(faceCentres, faceI)
    // {
    //   os << faceI << ","
    // 	 << faceCentres[faceI].x() << ","
    // 	 << faceCentres[faceI].y() << ","
    // 	 << faceCentres[faceI].z() << endl;
    // }

    // Return face area vector as surfaceVectorField
    const vectorField& faceAreaVectors = mesh.Sf().boundaryField()[patchI];
    // Write face area vectors to a file
    fileName outputFileSf("buildingPatchFaceAreaVectors");
    OFstream os2(outputFileSf);

    forAll(faceAreaVectors, faceI)
    {
        os2 << faceAreaVectors[faceI].x() << "," 
            << faceAreaVectors[faceI].y() << "," 
            << faceAreaVectors[faceI].z() << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
