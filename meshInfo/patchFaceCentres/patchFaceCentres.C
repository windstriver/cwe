#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    // Prepare options
    argList::addOption
    (
        "patchName",	
        "patch name where to extract patch face centres"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Create the output path directory
    fileName outputDir = mesh.time().path()/"constant/meshInfo";
    mkDir(outputDir);
    autoPtr<OFstream> outputFilePtr;

    /*************************************************************************
    Write inlet patch face centres
    *************************************************************************/

    if (args.optionFound("patchName"))
    {
        // Get the name of the patch to extract face centres
        const word patchName = args.option("patchName");

        // Find patch index given a patch name
        label patch = mesh.boundaryMesh().findPatchID(patchName);
        // Create face centres vectorField
        const vectorField& faceCentres = mesh.Cf().boundaryField()[patch];
        Info<< "Face centres of patch " << patchName << ":" << endl
            << faceCentres << endl;

        // Open the file from the output directory
        outputFilePtr.reset(new OFstream(outputDir+"/faceCentres"+patchName));

        forAll(faceCentres, faceI)
        {
            outputFilePtr() << "("
                            << faceCentres[faceI].x() << "    "
                            << faceCentres[faceI].y() << "    "
                            << faceCentres[faceI].z() << ")" << endl;
        }

        Info<< "Wrote face centres of patch " << patchName << " to constant/meshInfo/"
            << "faceCentres" << patchName << endl;

    }

    return 0;
}
