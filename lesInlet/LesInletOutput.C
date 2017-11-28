#include "LesInlet.H"
#include <iostream>
#include <vector>
#include <string>
#include "H5Cpp.h"

using namespace H5;

const std::string FILE_NAME("lesinlet.h5");
const std::string DATASET_NAME_U("U");
const std::string DATASET_NAME_V("V");
const std::string DATASET_NAME_W("W");
const std::string DATASET_NAME_GRID("GRID");
const std::string DATASET_NAME_TIME("TIME");
const int NX = 2;    // number of inlet mesh points
const int NY = 2000;    // number of time steps, including time 0
const int RANK = 2;
const double lastTime = 10;    // last time [sec]
const double dt = lastTime / (1.0*(NY-1));

int main()
{
    // Creat HDF5 file
    // H5File file(FILE_NAME, H5F_ACC_TRUNC);
    // hsize_t dims[2];
    // dims[0] = NX;
    // dims[1] = NY;
    // DataSpace dataspace(RANK, dims);
    // DataSet dataset_u = file.createDataSet(DATASET_NAME_U,
    //                                        PredType::IEEE_F64BE,
    //                                        dataspace);
    // DataSet dataset_v = file.createDataSet(DATASET_NAME_V,
    //                                        PredType::IEEE_F64BE,
    //                                        dataspace);
    // DataSet dataset_w = file.createDataSet(DATASET_NAME_W,
    //                                        PredType::IEEE_F64BE,
    //                                        dataspace);

    // dims[0] = NX;
    // dims[1] = 3;
    // DataSpace dataspace_grid(RANK, dims);
    // DataSet dataset_grid = file.createDataSet(DATASET_NAME_GRID,
    //                                           PredType::IEEE_F64BE,
    //                                           dataspace_grid);

    // hsize_t dims_time[1];
    // dims_time[0]= NY;
    // DataSpace dataspace_time(1, dims_time);
    // DataSet dataset_time = file.createDataSet(DATASET_NAME_TIME,
    //                                           PredType::IEEE_F64BE,
    //                                           dataspace_time);


    // Generate turbulence
    LesInlet inlet;
    // Store inlet mesh points
    vector<vector<double> > X(2, vector<double> (3));
    X[0][0] = 0;
    X[0][1] = 0;
    X[0][2] = 0.95;
    X[1][0] = 0;
    X[1][1] = 0;
    X[1][2] = 0.35;
    double grid[2][3];
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 3; ++j)
            grid[i][j] = X[i][j];
    // Store time steps vector
    vector<double> timeVec(NY, 0.0);
    double time[NY];
    for(int i = 0; i < NY; ++i)
    {
        timeVec[i] = i*dt;
        time[i] = timeVec[i];
    }

    // Turbulent velocity time history vector
    vector<vector<vector<double> > > UTh;
    UTh = inlet.Uturb(X, timeVec);
    double Ux[NX][NY];
    double Uy[NX][NY];
    double Uz[NX][NY];

    for(int i = 0; i < NX; ++i)
        for(int j = 0; j < NY; ++j)
        {
            Ux[i][j] = UTh[i][j][0];
            Uy[i][j] = UTh[i][j][1];
            Uz[i][j] = UTh[i][j][2];
        }


    // Write data to HDF5
    // Open an existing file and dataset.
    H5File file(FILE_NAME, H5F_ACC_RDWR);
    DataSet dataset = file.openDataSet(DATASET_NAME_U);
    dataset.write(Ux, PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet(DATASET_NAME_V);
    dataset.write(Uy, PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet(DATASET_NAME_W);
    dataset.write(Uz, PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet(DATASET_NAME_GRID);
    dataset.write(grid, PredType::NATIVE_DOUBLE);
    dataset = file.openDataSet(DATASET_NAME_TIME);
    dataset.write(time, PredType::NATIVE_DOUBLE);

    return 0;
}
