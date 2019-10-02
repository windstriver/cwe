/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

Application
    velocityFieldInit

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"
//#include "Kmesh.H"
//#include "turbGen.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Parameters
    scalar Ue = 6.355;
    scalar delta = 0.02835;
    scalar uTau = 0.2824;
    scalar nu = 1.4612e-5;
    scalar pi = Foam::constant::mathematical::pi;

    // Parameters for law of the wake (Coles, 1956)
    scalar kappa = 0.4;
    scalar B = 5.2;
    scalar PI = 0.3942;

    /*
    // turbGen from OpenFOAM
    scalar Ea = 2;
    scalar k0 = 1000;
    Kmesh K(mesh);
    Info<< "K field: " << K << endl;
    turbGen Ugen(K, Ea, k0);
    vectorField UPrime = Ugen.U();
    Info<< "Turbulence field: " << UPrime << endl;
    */

    // Mean velocity field based on law of the wake
    vectorField Umean(mesh.nCells(), vector(0., 0., 0.));
    vectorField UPrime(mesh.nCells(), vector(0., 0., 0.));

    Random RandGen(label(0));
    forAll(Umean, cellI)
    {
        scalar y = mesh.C()[cellI].y();
        if (y>delta) {y=delta;}
        scalar yPlus = y*uTau/nu;
        scalar eta = y/delta;
        Umean[cellI].x() = uTau * (1.0/kappa*Foam::log(yPlus)+B+2*PI/kappa*Foam::sin(pi/2*eta)*Foam::sin(pi/2*eta));
        UPrime[cellI] = (RandGen.sample01<vector>()-vector(0.5, 0.5, 0.5))*0.2*Ue;
    }

    vectorField U = Umean + UPrime;
    // Write generated velocity field to a file in 0 dir
    Info<< "timePath: " << runTime.timePath() << endl;

    fileName outputDir = runTime.timePath();
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr.reset(new OFstream(outputDir/"Uinit"));
    writeEntry(outputFilePtr(), "", U);
    return 0;
}

// ************************************************************************* //
