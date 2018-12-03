#include "fvCFD.H"
#include "H5Cpp.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    // Determine the number of inlet patch face centers
    Foam::label patchI = mesh.boundaryMesh().findPatchID("inlet");
    const int patchFaceNum = mesh.boundaryMesh()[patchI].faceCentres().size();
    Foam::Info <<
        "Total number of faces at the inlet boundary patch: " << patchFaceNum << Foam::endl;
    // Prepare the inlet patch face centers data
    double patchFaceCenter[patchFaceNum][4];
    int i = 0;
    forAll(mesh.boundaryMesh()[patchI].faceCentres(), faceI)
    {
        patchFaceCenter[i][0] = faceI;
        patchFaceCenter[i][1] = mesh.boundaryMesh()[patchI].faceCentres()[faceI].x();
        patchFaceCenter[i][2] = mesh.boundaryMesh()[patchI].faceCentres()[faceI].y();
        patchFaceCenter[i][3] = mesh.boundaryMesh()[patchI].faceCentres()[faceI].z();
        ++i;
    }

    // Determine the time steps
    Foam::IOdictionary controlDict
    (
        Foam::IOobject
        (
            "controlDict",
            runTime.system(),
            mesh,
            Foam::IOobject::MUST_READ_IF_MODIFIED,
            Foam::IOobject::NO_WRITE
        )
    );

    double startTime = controlDict.lookup("startTime")[0].scalarToken();
    double endTime = controlDict.lookup("endTime")[0].scalarToken();
    double deltaT = controlDict.lookup("deltaT")[0].scalarToken();
    Foam::Info << "startTime: " << startTime << Foam::endl;
    Foam::Info << "endTime: " << endTime << Foam::endl;
    Foam::Info << "deltaT: " << deltaT << Foam::endl;

    // number of inlet mesh points
    const int NX = patchFaceNum;
    // number of time steps, including time 0
    const int NT = int((endTime-startTime)/deltaT+1);

    // Store time steps vector
    double time[NT];
    for(int i = 0; i < NT; ++i)
    {
        time[i] = startTime + i*deltaT;
    }


    // Parameters for HDF5 databases
    const std::string FILE_NAME("../inflowTurbCDRFG.h5");
    const std::string DATASET_NAME_UMEAN("UMEAN");
    const std::string DATASET_NAME_U("U");
    const std::string DATASET_NAME_V("V");
    const std::string DATASET_NAME_W("W");
    const std::string DATASET_NAME_GRID("GRID");
    const std::string DATASET_NAME_TIME("TIME");

    // Creat HDF5 file
    H5::H5File file(FILE_NAME, H5F_ACC_TRUNC);
    // Uturb
    hsize_t dims[2];
    dims[0] = NX;
    dims[1] = NT;
    H5::DataSpace dataspace(2, dims);
    H5::DataSet dataset_u = file.createDataSet(DATASET_NAME_U,
                                               H5::PredType::IEEE_F64BE,
                                               dataspace);
    H5::DataSet dataset_v = file.createDataSet(DATASET_NAME_V,
                                               H5::PredType::IEEE_F64BE,
                                               dataspace);
    H5::DataSet dataset_w = file.createDataSet(DATASET_NAME_W,
                                               H5::PredType::IEEE_F64BE,
                                               dataspace);
    // Umean
    hsize_t dims_umean[1];
    dims_umean[0] = NX;
    H5::DataSpace dataspace_umean(1, dims_umean);
    H5::DataSet dataset_umean = file.createDataSet(DATASET_NAME_UMEAN,
                                                   H5::PredType::IEEE_F64BE,
                                                   dataspace_umean);
    // Grid
    dims[0] = NX;
    dims[1] = 4;
    H5::DataSpace dataspace_grid(2, dims);
    H5::DataSet dataset_grid = file.createDataSet(DATASET_NAME_GRID,
                                                  H5::PredType::IEEE_F64BE,
                                                  dataspace_grid);
    // Time
    hsize_t dims_time[1];
    dims_time[0]= NT;
    H5::DataSpace dataspace_time(1, dims_time);
    H5::DataSet dataset_time = file.createDataSet(DATASET_NAME_TIME,
                                                  H5::PredType::IEEE_F64BE,
                                                  dataspace_time);

    // Write mesh and time to HDF5
    H5::H5File file2(FILE_NAME, H5F_ACC_RDWR);
    H5::DataSet dataset2 = file2.openDataSet(DATASET_NAME_GRID);
    dataset2.write(patchFaceCenter, H5::PredType::NATIVE_DOUBLE);
    dataset2 = file2.openDataSet(DATASET_NAME_TIME);
    dataset2.write(time, H5::PredType::NATIVE_DOUBLE);

    return 0;
}


// ************************************************************************* //
