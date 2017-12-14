#ifndef LesInlet_H
#define LesInlet_H

#include <cmath>
#include <random>
#include <vector>

#include "WindProfile.H"

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
    std::vector<std::vector<std::vector<double> > >
    UturbFreqSeg(double freq,
                 std::vector<std::vector<double> > xMat,
                 std::vector<double> tt);

    // Veclotiy for whole frequency domain
    // Input:
    //    xMat: matrix to store inlet point c.s. nPts x 3 std::vector
    //    tt:   std::vector to store time step std::vector
    // Output:
    //    Matrix to store turbulent velocity at inlet point and time step
    std::vector<std::vector<std::vector<double> > >
    Uturb(std::vector<std::vector<double> > xMat,
          std::vector<double> tt);

};

#endif
