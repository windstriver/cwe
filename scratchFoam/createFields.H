Info<< "Reading field T\n" << endl;

    // Create scalar field T
    volScalarField T
    (
        // Create object for input/output operations
        IOobject
        (
            // Name of the dictionary file to read/write
            "T",
            // runtime directory
            runTime.timeName(),
            // Object registry
            mesh,
            // Read the directory in the runtime directory
            // (MUST_READ)
            IOobject::MUST_READ,
            // Write the value in the runtime directory
            // (AUTO_WRITE)
            // If you do not want to write the value,
            // use the option NO_WRITE.
            IOobject::AUTO_WRITE
        ),
        mesh    // Link object to mesh
    );

    // We add the initialization of T after its declaration
    forAll (T, i)
    {
        // Access cell center coordinates
        const scalar x = mesh.C()[i][0];
        const scalar y = mesh.C()[i][1];
        const scalar z = mesh.C()[i][2];

        if ( 0.3 < x && x < 0.7 )
        {
            T[i] = 1.0;
        }
    }
    T.write();    // Write field T. As the file createFields.H is
                  // outside of the time loop the value is saved
                  // in the time directory 0


Info<< "Reading field U\n" << endl;

    // Create vector field U
    volVectorField U
    (
        IOobject
        (
            // Name of the dictionary file to read/write
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


Info<< "Reading transportProperties\n" << endl;

    // Create object transportProperties used to read data
    IOdictionary transportProperties
    (
        IOobject
        (
            // Name of the input dictionary
            "transportProperties",
            // Location of the input dictionary
            // in this case is located in the directory "constant"
            runTime.constant(),
            mesh,
            // Re-read data if it is modified
            IOobject::MUST_READ_IF_MODIFIED,
            // Do not write anything in the dictionary
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading diffusivity DT\n" << endl;

    // Create scalar DT (diffusion coefficient)
    dimensionedScalar DT
    (
        // Access value of DT in the object transportProperties
        transportProperties.lookup("DT")
    );


    // Create and initialize the relative face flux field phi
    #include "createPhi.H"

