#ifndef LesInlet_H
#define LesInlet_H

#include <cmath>
#include <random>
#include <vector>

#include "WindProfile.H"

using namespace std;

class LesInlet: public WindProfile
{
    // Private data

    // Number of frequency segments
    int nm;
    // Number of random frequencies in one segment
    int nf;
    // Maximum frequency
    double fmax;
    double fmin;
    double df;

public:

    // Constructors
    LesInlet();

    vector<double> Um(double freq, vector<double> X, double t);

};

#endif
