#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Create the output path directory
    fileName outputDir = mesh.time().path()/"constant/meshInfo";
    mkDir(outputDir);
    autoPtr<OFstream> outputFilePtr;

    /*************************************************************************
    Write cell centres in a box region
    *************************************************************************/

    // Find patch index given a patch name
    // Create face centres vectorField
    const vectorField& cellCentres = mesh.Cf().internalField();

    // Open the file from the output directory
    outputFilePtr.reset(new OFstream(outputDir+"/cellCentresMidPlane"));

    forAll(cellCentres, cellI)
    {
        if (cellCentres[cellI].x() > (0.15-1.0e-3) && 
            cellCentres[cellI].x() < (0.15+1.0e-3) )
        {
            outputFilePtr() << "("
                            << cellCentres[cellI].x() << "    "
                            << cellCentres[cellI].y() << "    "
                            << cellCentres[cellI].z() << ")" << endl;
        }
    }

    Info<< "Wrote cell centres near mid-plane to constant/meshInfo/cellCentresMidPlane" << endl;

    return 0;
}
