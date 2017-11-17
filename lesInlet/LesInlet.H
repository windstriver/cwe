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
    vector<double> UturbFreqSeg(double freq, vector<double> X, double t);

    // Veclotiy for whole frequency domain
    vector<double> Uturb(vector<double> X, double t);

};

#endif
