#ifndef LesInlet_H
#define LesInlet_H

#include <cmath>
#include <random>
#include <vector>

#include "WindProfile.H"

using namespace std;

class LesInlet: WindProfile
{
    // Private data

    // Number of frequency segments
    int nm;
    // Number of random frequencies in one segment
    int nf;
    // Maximum frequency
    double fmax;
    // Minimum frequency
    double fmin;
    // Delta frequency
    double df;

public:

    // Constructors
    LesInlet();

    // Access Maximum frequency
    double getFmax();

    // Velocity for frequency single frequency
    vector<vector<vector<double> > >
    UturbFreqSeg(double freq,
                 vector<vector<double> > xMat,
                 vector<double> tt);

    // Veclotiy for whole frequency domain
    // Input:
    //    xMat: matrix to store inlet point c.s. nPts x 3 vector
    //    tt:   vector to store time step vector
    // Output:
    //    Matrix to store turbulent velocity at inlet point and time step
    vector<vector<vector<double> > >
    Uturb(vector<vector<double> > xMat,
          vector<double> tt);

};

#endif
