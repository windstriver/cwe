#include <fstream>
#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    Foam::label patchI = mesh.boundaryMesh().findPatchID("inlet");

    std::ofstream file;
    file.open("inletCenter.csv");
    forAll(mesh.boundaryMesh()[patchI].faceCentres(), faceI)
    {
        Foam::scalar x = mesh.boundaryMesh()[patchI].faceCentres()[faceI].x();
        Foam::scalar y = mesh.boundaryMesh()[patchI].faceCentres()[faceI].y();
        Foam::scalar z = mesh.boundaryMesh()[patchI].faceCentres()[faceI].z();
        file << faceI << "," << x << "," << y << "," << z << "," << std::endl;
    }
    file.close();

    return 0;
}
