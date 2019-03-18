#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    // Building height
    scalar H = 0.364;
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Create cell centres vectorField
    const vectorField& cellCentres = mesh.C().internalField();

    // Create the output path directory
    fileName outputDir = mesh.time().path()/"constant/meshInfo";
    mkDir(outputDir);

    /*************************************************************************
    Write cell centres
    *************************************************************************/
    // Write cell centres to file
    // File pointer to direct the output to
    autoPtr<OFstream> outputFilePtr;
    autoPtr<OFstream> outputFilePtr2;
    // Open the file from the output directory
    outputFilePtr.reset(new OFstream(outputDir/"cellCentres.csv"));
    outputFilePtr() << "cellNo , probeNo , x , y , z" << endl;
    outputFilePtr2.reset(new OFstream(outputDir/"cellCentresOF"));

    int i = 0;
    scalar x, y, z;
    forAll(cellCentres, cellI)
    {
	z = cellCentres[cellI].z();
	y = cellCentres[cellI].y();
	if ((z < H) && (y < H) && (y > -1.0*H))
	{
	    x = cellCentres[cellI].x();

	    outputFilePtr() << cellI << ","
			    << i     << ","
			    << x     << ","
			    << y     << ","
			    << z     << endl;
	    
	    outputFilePtr2() << "("
			     << x     << "    "
			     << y     << "    "
			     << z     << ")" << endl;
	    i++;
	}
    }
    
    /*************************************************************************
    Write inlet patch face centres
    *************************************************************************/

    // Find patch index given a patch name
    label patch = mesh.boundaryMesh().findPatchID("inlet");
    // Create face centres vectorField
    const vectorField& faceCentres = mesh.Cf().boundaryField()[patch];
    // Open the file from the output directory
    outputFilePtr.reset(new OFstream(outputDir/"faceCentresInlet.csv"));
    outputFilePtr2.reset(new OFstream(outputDir/"faceCentresInletOF"));

    forAll(faceCentres, faceI)
    {
	x = faceCentres[faceI].x();
	y = faceCentres[faceI].y();
	z = faceCentres[faceI].z();

	outputFilePtr()	<< faceI << ","
			<< x     << ","
			<< y     << ","
			<< z     << endl;
	    
	outputFilePtr2() << "("
			 << x     << "    "
			 << y     << "    "
			 << z     << ")" << endl;
    }

    /*************************************************************************
    Write inlet patch face area vectors
    *************************************************************************/

    const vectorField& faceAreaVectors = mesh.Sf().boundaryField()[patch];
    // Write face area vectors to a file
    outputFilePtr.reset(new OFstream(outputDir/"faceAreaVectorsInlet.csv"));
    forAll(faceAreaVectors, faceI)
    {
        outputFilePtr() << faceAreaVectors[faceI].x() << "," 
                        << faceAreaVectors[faceI].y() << "," 
                        << faceAreaVectors[faceI].z() << endl;
    }

    /*************************************************************************
    Write ground patch face centres
    *************************************************************************/

    // Find patch index given a patch name
    patch = mesh.boundaryMesh().findPatchID("ground");
    // Create face centres vectorField
    const vectorField& faceCentres2 = mesh.Cf().boundaryField()[patch];
    // Open the file from the output directory
    outputFilePtr.reset(new OFstream(outputDir/"faceCentresGround.csv"));
    outputFilePtr() << "cellNo , probeNo , x , y" << endl;
    outputFilePtr2.reset(new OFstream(outputDir/"faceCentresGroundOF"));

    i = 0;
    forAll(faceCentres2, faceI)
    {
	z = faceCentres2[faceI].z();
	y = faceCentres2[faceI].y();
	if ((y < H) && (y > -1.0*H))
	{
	    x = faceCentres2[faceI].x();

	    outputFilePtr() << faceI << ","
			    << i     << ","
			    << x     << ","
			    << y     << ","
			    << z     << endl;
	    
	    outputFilePtr2() << "("
			     << x     << "    "
			     << y     << "    "
			     << z     << ")" << endl;
	    i++;
	}
    }
    
    return 0;

}
