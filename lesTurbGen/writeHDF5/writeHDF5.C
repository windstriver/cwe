#include <vector>
#include <iostream>
// #include "fvCFD.H"
#include "H5Cpp.h"
#include "LesInlet.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Parameters for HDF5 databases
    const std::string FILE_NAME("../lesinlet.h5");
    const std::string DATASET_NAME_U("U");
    const std::string DATASET_NAME_V("V");
    const std::string DATASET_NAME_W("W");
    const std::string DATASET_NAME_GRID("GRID");
    const std::string DATASET_NAME_TIME("TIME");

    // Get grid size
    H5::H5File file(FILE_NAME, H5F_ACC_RDWR);
    H5::DataSet dataset = file.openDataSet(DATASET_NAME_GRID);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims_grid[2];
    int ndims = dataspace.getSimpleExtentDims(dims_grid, NULL);
    const int NX = dims_grid[0];
    // Define the memory dataspace
    hsize_t dims_mem_grid[2];
    dims_mem_grid[0] = NX;
    dims_mem_grid[1] = 4;
    H5::DataSpace mem_grid_space(2, dims_mem_grid);
    // Read grid c.s.
    double grid[NX][4];
    dataset.read(grid, H5::PredType::NATIVE_DOUBLE, mem_grid_space, dataspace);
    // for (int i = 0; i < NX; ++i)
    // {
    //     for (int j = 0; j < 4; ++j)
    //     {
    //         std::cout << grid[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }


    // Get time steps
    dataset = file.openDataSet(DATASET_NAME_TIME);
    dataspace = dataset.getSpace();
    hsize_t dims_time[1];
    ndims = dataspace.getSimpleExtentDims( dims_time, NULL);
    const int NT = dims_time[0];
    // std::cout << NT << std::endl;
    // Define the memory dataspace
    hsize_t dims_mem_time[1];
    dims_mem_time[0] = NT;
    H5::DataSpace mem_time_space(1, dims_mem_time);
    // Read time steps
    double time[NT];
    dataset.read(time, H5::PredType::NATIVE_DOUBLE, mem_time_space, dataspace);
    // for (int i = 0; i < NT; ++i)
    // {
    //     std::cout << time[i] << std::endl;
    // }

    // Convert mesh grids and time steps to vector
    std::vector<std::vector<double> > XVec(NX, std::vector<double> (3));
    for (int i = 0; i < NX; ++i)
    {
        XVec[i][0] = grid[i][3];
        XVec[i][1] = grid[i][1];
        XVec[i][2] = grid[i][2];
    }
    std::vector<double> timeVec(NT, 0.0);
    for (int i = 0; i < NT; ++i)
    {
        timeVec[i] = time[i];
    }

    LesInlet inlet;
    std::vector<std::vector<std::vector<double> > > UTh;
    UTh = inlet.Uturb(XVec, timeVec);
    double Ux[NX][NT];
    double Uy[NX][NT];
    double Uz[NX][NT];

    for(int i = 0; i < NX; ++i)
    {
        for(int j = 0; j < NT; ++j)
        {
            Ux[i][j] = UTh[i][j][0];
            Uy[i][j] = UTh[i][j][1];
            Uz[i][j] = UTh[i][j][2];
            // std::cout << UTh[i][j][0] << std::endl;
        }
    }

    // Write data to HDF5
    // Open an existing file and dataset.
    dataset = file.openDataSet(DATASET_NAME_U);
    dataset.write(Ux, H5::PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet(DATASET_NAME_V);
    dataset.write(Uy, H5::PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet(DATASET_NAME_W);
    dataset.write(Uz, H5::PredType::NATIVE_DOUBLE);

    return 0;
}


// ************************************************************************* //
