#include "LesInlet.H"
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

int main()
{
    WindProfile wp;
    LesInlet inlet;

    // Point at 0.95 m height
    vector<vector<double> > X(2, vector<double> (3));
    X[0][0] = 0;
    X[0][1] = 0;
    X[0][2] = 0.95;
    X[1][0] = 0;
    X[1][1] = 0;
    X[1][2] = 0.35;

    // Test of von Karman spectrum
    double fmax = inlet.getFmax();
    double fmin = 1;
    int nf = 100;
    double df = (fmax-fmin) / (nf-1);
    // Frequency vector
    vector<double> freqVec(nf, 0.0);
    vector<double> SuVec(nf, 0.0);
    for (int i = 0; i < nf; ++i)
    {
        freqVec[i] = fmin + i*df;
        SuVec[i] = wp.vonKarmanSu(X[0][2], freqVec[i]);
        // cout << "Freq: " << freqVec[i] << "    " << "Su: " << SuVec[i] << endl;
    }

    ofstream output;
    // Creat a file
    output.open("Su.csv");
    for (int i = 0; i < nf; ++i)
    {
        output << freqVec[i] << "," << SuVec[i] << endl;
    }
    output.close();


    // Number of time steps
    int nt = 2000;
    // Time step
    double dt = 1.0/fmax/2/2.5;
    // Time vector
    vector<double> timeVec(nt, 0.0);
    for(int i = 0; i < nt; ++i)
    {
        timeVec[i] = i*dt;
    }

    // Turbulent velocity time history vector
    vector<vector<vector<double> > > UTh;
    UTh = inlet.Uturb(X, timeVec);
    vector<double> UxTh1(timeVec.size());
    // vector<double> UxTh2;

    for (unsigned int ts = 0; ts < timeVec.size(); ++ts)
    {
        UxTh1[ts] = UTh[0][ts][0];
        // cout << UxTh1[ts] << endl;
    }

    // Calculate the mean and standard deviation of Ux time history
    double meanUx = 0;
    double stdUx = 0;
    for (int i = 0; i < nt; ++i)
    {
        meanUx += UxTh1[i];
    }
    meanUx /= nt;

    for (int i = 0; i < nt; ++i)
    {
        stdUx += (UxTh1[i] - meanUx) * (UxTh1[i] - meanUx);
    }
    stdUx /= (nt-1);

    cout << "Mean of Ux: " << meanUx << endl
         << "Std  of Ux: " << stdUx << endl;

    // Creat a file
    output.open("UxTh.csv");
    for (int i = 0; i < nt; ++i)
    {
        output << timeVec[i] << "," << UxTh1[i] << endl;
    }
    output.close();

    return 0;
}
