#ifndef WindProfile_H
#define WindProfile_H

class WindProfile
{
    // Private data

    // Mean Velocity
    // Reference height for the mean velocity
    double h0u;
    // Mean velocity at h0u
    double Uh;
    // Power low exponent of the mean velocity
    double alphau;

    // Turbulence Intensity
    // Reference height for the turbulenct intensity
    double h0I;
    // Longitudinal turbulent intensity at h0I
    double Iuh;
    // Transverse turbulent intensity at h0I
    double Ivh;
    // Vertical turbulent intensity at h0I
    double Iwh;
    // Power law exponent of the longitudinal turbulent intensity
    double dIu;
    // Power law exponeent of the transverse turbulent intensity
    double dIv;
    // Power law exponent of the vertical turbulent intensity
    double dIw;

    // Turbulence length scale
    // Reference height for the length scale
    double h0L;
    // Longitudinal length scale at h0L
    double Luh;
    // Transverse length scale at h0L
    double Lvh;
    // Vertical length scale at h0L
    double Lwh;
    // Power law exponent of the longitudinal length scale
    double dLu;
    // Power law exponent of the transverse length scale
    double dLv;
    // Power law exponent of the vertical length scale
    double dLw;

    // Coherency
    // Coherency decay constants in x, y and z directions
    double Cx;
    double Cy;
    double Cz;
    // Characteristic length to maintain the coherency
    double dGamma;

public:

    // Constructors
    WindProfile();

    // Member Functions
    // Calculate mean velocity
    double meanVelocity(double z);

    // Calculate  turbulent intensity
    double turbIntensityU(double z);
    double turbIntensityV(double z);
    double turbIntensityW(double z);

    // Calculate turbulence length scale
    double turbLengthScaleU(double z);
    double turbLengthScaleV(double z);
    double turbLengthScaleW(double z);

    // von Karman turbulent spectra
    double vonKarmanSu(double z, double freq);
    double vonKarmanSv(double z, double freq);
    double vonKarmanSw(double z, double freq);

    // Calculate tuning factor
    double tuningFactor(double z);

    // Acess the coherency decay constant
    double getCx();
    double getCy();
    double getCz();

};

#endif
