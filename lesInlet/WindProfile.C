#include "WindProfile.H"
#include <cmath>

WindProfile::WindProfile()
{
    // Mean velocity
    h0u = 0.3644;    // [m]
    Uh = 10.0;      // [m/s]
    alphau = 0.3264;

    // Turbulent intensity
    h0I = 0.3364;    // [m]
    Iuh = 0.2084;
    Ivh = 0.1815;
    Iwh = 0.1523;
    dIu = 0.1914;
    dIv = 0.1228;
    dIw = 0.0048;

    // Turbulence length scale
    h0L = 0.254;   // [m]
    Luh = 0.302;    // [m]
    Lvh = 0.0815;   // [m]
    Lwh = 0.0326;   // [m]
    dLu = 0.473;
    dLv = 0.8813;
    dLw = 1.5390;

    // Coherency
    Cx = 10;
    Cy = 10;
    Cz = 10;
    dGamma = 0.3;

}


double WindProfile::meanVelocity(double z)
{
    return Uh * pow((z/h0u), alphau);
}

double WindProfile::turbIntensityU(double z)
{
    return Iuh * pow((z/h0I), -1.0*dIu);
}

double WindProfile::turbIntensityV(double z)
{
    return Ivh * pow((z/h0I), -1.0*dIv);
}

double WindProfile::turbIntensityW(double z)
{
    return Iwh * pow((z/h0I), -1.0*dIw);
}

double WindProfile::turbLengthScaleU(double z)
{
    return Luh * pow((z/h0L), dLu);
}

double WindProfile::turbLengthScaleV(double z)
{
    return Lvh * pow((z/h0L), dLv);
}

double WindProfile::turbLengthScaleW(double z)
{
    return Lwh * pow((z/h0L), dLw);
}

double WindProfile::vonKarmanSu(double z, double freq)
{
    double Uav = meanVelocity(z);
    double Iu = turbIntensityU(z);
    double Lu = turbLengthScaleU(z);
    return 4*pow(Iu*Uav,2)*(Lu/Uav) /                \
        pow((1+70.8*pow(freq*Lu/Uav,2)), 5.0/6.0);
}

double WindProfile::vonKarmanSv(double z, double freq)
{
    double Uav = meanVelocity(z);
    double Iv = turbIntensityV(z);
    double Lv = turbLengthScaleV(z);
    return 4*pow(Iv*Uav,2)*(Lv/Uav) * (1+188.4*pow(2*freq*Lv/Uav,2)) /  \
        pow((1+70.8*pow(2*freq*Lv/Uav,2)), 11.0/6.0);
}

double WindProfile::vonKarmanSw(double z, double freq)
{
    double Uav = meanVelocity(z);
    double Iw = turbIntensityV(z);
    double Lw = turbLengthScaleV(z);
    return 4*pow(Iw*Uav,2)*(Lw/Uav) * (1+188.4*pow(2*freq*Lw/Uav,2)) /  \
        pow((1+70.8*pow(2*freq*Lw/Uav,2)), 11.0/6.0);
}

double WindProfile::gamma(double z)
{
    double Lu = turbLengthScaleU(z);
    double beta = Cx * dGamma / Lu;
    if (beta<6.0)
        return 3.7 * pow(beta,-0.3);
    else
        return 2.1;
}
