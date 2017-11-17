#include "LesInlet.H"
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

int main()
{
    LesInlet wp;

    double fmax = wp.getFmax();

    // Point at 0.95 m height
    vector<double> X(3);
    X[0] = 0;
    X[1] = 0;
    X[2] = 0.95;


    // Number of time steps
    int nt = 2000;
    // Time step
    double dt = 1.0/fmax/2/2.5;
    // Time vector
    vector<double> timeVec(nt, 0.0);

    // Turbulent velocity time history vector
    vector<double> UxTh(nt, 0.0);
    vector<double> UyTh(nt, 0.0);
    vector<double> UzTh(nt, 0.0);

    for(int i = 0; i < nt; ++i)
    {
        timeVec[i] = i*dt;
        vector<double> U = wp.Uturb(X, timeVec[i]);
        UxTh[i] = U[0];
        UyTh[i] = U[1];
        UzTh[i] = U[2];
        cout << "Time step: " << timeVec[i] << "    Ux: " << UxTh[i] << endl;
    }

    // Calculate the mean and standard deviation of Ux time history
    double meanUx = 0;
    double stdUx = 0;
    for (int i = 0; i < nt; ++i)
    {
        meanUx += UxTh[i];
    }
    meanUx /= nt;

    for (int i = 0; i < nt; ++i)
    {
        stdUx += (UxTh[i] - meanUx) * (UxTh[i] - meanUx);
    }
    stdUx /= (nt-1);

    cout << "Mean of Ux: " << meanUx << endl
         << "Std  of Ux: " << stdUx << endl;

    ofstream output;
    // Creat a file
    output.open("UxTh.csv");
    for (int i = 0; i < nt; ++i)
    {
        output << timeVec[i] << "," << UxTh[i] << endl;
    }
    output.close();


    return 0;
}
