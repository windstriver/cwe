volVectorField gradT(fvc::grad(T));

volVectorField gradT_vector
(
	IOobject
	(
		"gradT",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	gradT
);

volScalarField gradTx
(
	IOobject
	(
		"gradTx",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	gradT.component(vector::X)
	//gradT.component(0)
);

volScalarField gradTy
(
	IOobject
	(
		"gradTy",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	gradT.component(vector::Y)
);

volScalarField gradTz
(
	IOobject
	(
		"gradTz",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
	gradT.component(vector::Z)
);


volScalarField divGradT
(
	IOobject
	(
		"divGradT",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::div(gradT)
);


dimensionedScalar myDT
(
	"myDt",
	//dimensionSet(0, 2, 0, 0, 0, 0, 0),
    dimless,
	1
);

volScalarField my_laplacianT
(
	IOobject
	(
		"my_laplacianT",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::laplacian(myDT, T)
);


gradT_vector.write();
gradTx.write();
gradTy.write();
gradTz.write();
divGradT.write();
my_laplacianT.write();
