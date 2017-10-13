/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    scratchFoam

Description

\*---------------------------------------------------------------------------*/

// Important header
// add all the class declarations needed to access mesh, fields,
// tensor algebra, fvm/fvc operations, time, parallel communication,
// linear algebra, and so on.
#include "fvCFD.H"

// Solution control using PISO class
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Set directory structure
    #include "setRootCase.H"

    // Create time (object runtime)
    #include "createTime.H"

    // Create mesh (object mesh)
    #include "createMesh.H"

    // Initialize fields (need to create such a file)
    #include "createFields.H"

    // Calculates and outputs the Courant number
    #include "CourantNo.H"

    // Declare and initialize the cumulative continuity error
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Assign PISO controls to object mesh
    pisoControl piso(mesh);

    // Output some information
    Info<< "\nStarting time loop\n" << endl;

    // Time loop
    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        while (piso.correct())    // PISO options (current loop)
        {
            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix TEqn    // Create object Teqn.
                                       // fvScalarMatrix is a scalar instance
                                       // of fvMatrix
                (
                    fvm::ddt(T)
                    + fvm::div(phi, T)
                    - fvm::laplacian(DT, T)

                    // Model equation (convection-diffusion)
                    // need to create the scalar field T, vector field U
                    // (used in phi or face flux), and the constant DT.
                    // We will declare these variables in the
                    // createFields.H header file.
                    // In the dictionary fvSchemes, you will need to
                    // define how to compute the differential operators,
                    // that is,
                    //     ddt(T)
                    //     div(phi, T)
                    //     laplacian(DT, T)
                );

                TEqn.solve();
                // Solve TEqn
                // At this point, the object TEqn holds the solution

            }
        }

        // Computes continuity errors
        #include "continuityErrs.H"

        // Add this header file to write extra fields
        #include "write.H"

        // Write the solution in the runtime folder
        // it will write the data requested in the file createFields.H
        runTime.write();
    }
    // At this point, we are outside of the time loop

    // Write CPU time at the end the time loop.
    // If you want to  compute the CPU time of each iteration,
    // add the same statement inside the time loop.
    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
    // End of the program (exit status)
    // If everything went fine, the program should return 0.
    // To now the return value, type in the terminal,
    // $> echo $?
}


// ************************************************************************* //
