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
        "patch name where to extract patch face area vectors"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Create the output path directory
    fileName outputDir = mesh.time().path()/"constant/meshInfo";
    mkDir(outputDir);
    autoPtr<OFstream> outputFilePtr;

    /*************************************************************************
    Write patch face area vectors
    *************************************************************************/

    if (args.optionFound("patchName"))
    {
        // Get the name of the patch to extract face centres
        const word patchName = args.option("patchName");

        // Find patch index given a patch name
        label patch = mesh.boundaryMesh().findPatchID(patchName);
        // Create face centres vectorField
        const vectorField& faceAreaVectors = mesh.Sf().boundaryField()[patch];
        Info<< "Face area vectors of patch " << patchName << ":" << endl
            << faceAreaVectors << endl;

        // Open the file from the output directory
        outputFilePtr.reset(new OFstream(outputDir+"/faceAreaVectors"+patchName));

        forAll(faceAreaVectors, faceI)
        {
            outputFilePtr() << "("
                            << faceAreaVectors[faceI].x() << "    "
                            << faceAreaVectors[faceI].y() << "    "
                            << faceAreaVectors[faceI].z() << ")" << endl;
        }

        Info<< "Wrote face centres of patch " << patchName << " to constant/meshInfo/"
            << "faceAreaVectors" << patchName << endl;

    }

    return 0;
}
