#include "WindProfile.H"
#include <cmath>

WindProfile::WindProfile()
{
    // Mean velocity
    h0u = 0.5;    // [m]
    Uh = 11.11;      // [m/s]
    alphau = 0.25;

    // Turbulent intensity
    h0I = 0.5;    // [m]
    Iuh = 0.116;
    Ivh = 0.087;
    Iwh = 0.058;
    dIu = 0;
    dIv = 0;
    dIw = 0;

    // Turbulence length scale
    h0L = 0.5;   // [m]
    Luh = 0.8 * 0.5;    // [m]
    Lvh = 0.8 * 0.5;   // [m]
    Lwh = 0.8 * 0.5;   // [m]
    dLu = 0;
    dLv = 0;
    dLw = 0;

    // Coherency
    Cx = 10;
    Cy = 10;
    Cz = 10;
    dGamma = 0.3;
}


double WindProfile::meanVelocity(double z)
{
    return Uh * std::pow((z/h0u), alphau);
}

double WindProfile::turbIntensityU(double z)
{
    return Iuh * std::pow((z/h0I), -1.0*dIu);
}

double WindProfile::turbIntensityV(double z)
{
    return Ivh * std::pow((z/h0I), -1.0*dIv);
}

double WindProfile::turbIntensityW(double z)
{
    return Iwh * std::pow((z/h0I), -1.0*dIw);
}

double WindProfile::turbLengthScaleU(double z)
{
    return Luh * std::pow((z/h0L), dLu);
}

double WindProfile::turbLengthScaleV(double z)
{
    return Lvh * std::pow((z/h0L), dLv);
}

double WindProfile::turbLengthScaleW(double z)
{
    return Lwh * std::pow((z/h0L), dLw);
}

double WindProfile::vonKarmanSu(double z, double freq)
{
    double Uav = meanVelocity(z);
    double Iu = turbIntensityU(z);
    double Lu = turbLengthScaleU(z);
    return 4*std::pow(Iu*Uav,2)*(Lu/Uav) /                \
        std::pow((1+70.8*std::pow(freq*Lu/Uav,2)), 5.0/6.0);
}

double WindProfile::vonKarmanSv(double z, double freq)
{
    double Uav = meanVelocity(z);
    double Iv = turbIntensityV(z);
    double Lv = turbLengthScaleV(z);
    return 4*std::pow(Iv*Uav,2)*(Lv/Uav) * (1+188.4*std::pow(2*freq*Lv/Uav,2)) / \
        std::pow((1+70.8*std::pow(2*freq*Lv/Uav,2)), 11.0/6.0);
}

double WindProfile::vonKarmanSw(double z, double freq)
{
    double Uav = meanVelocity(z);
    double Iw = turbIntensityV(z);
    double Lw = turbLengthScaleV(z);
    return 4*std::pow(Iw*Uav,2)*(Lw/Uav) * (1+188.4*std::pow(2*freq*Lw/Uav,2)) / \
        std::pow((1+70.8*std::pow(2*freq*Lw/Uav,2)), 11.0/6.0);
}

double WindProfile::tuningFactor(double z)
{
    double Lu = turbLengthScaleU(z);
    double beta = Cx * dGamma / Lu;
    if (beta<6.0)
        return 3.7 * std::pow(beta,-0.3);
    else
        return 2.1;
}

double WindProfile::getCx()
{
    return Cx;
}

double WindProfile::getCy()
{
    return Cy;
}

double WindProfile::getCz()
{
    return Cz;
}
