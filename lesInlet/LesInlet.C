#include "LesInlet.H"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

LesInlet::LesInlet()
{
    // Number of frequency segments
    nm = 50;
    // Number of random frequencies in one segment
    nf = 100;
    // Maximum frequency
    fmax = 100;
    // Minimum frequency
    fmin = fmax / nm / 2.0;
    // Frequency increment
    df = (fmax - fmin) / (nm - 1);
}


vector<double> LesInlet::Um(double freq, vector<double> X, double t)
{
    random_device rd;
    mt19937 mt(rd());
    normal_distribution<double> dist_r(0.0, 1.0);
    normal_distribution<double> dist_n(freq, df);

    vector<double> rx(nf);
    vector<double> ry(nf);
    vector<double> rz(nf);

    for(int i = 0; i < nf; ++i)
    {
        rx[i] = dist_r(mt);
        ry[i] = dist_r(mt);
        rz[i] = dist_r(mt);
    }

    vector<double> fmn(nf);
    for(int i = 0; i < nf; ++i)
    {
        fmn[i] = dist_n(mt);
    }

    vector<double> px(nf);
    vector<double> py(nf);
    vector<double> pz(nf);

    vector<double> qx(nf);
    vector<double> qy(nf);
    vector<double> qz(nf);

    for(int i = 0; i < nf; ++i)
    {
        // Calculate px, qx
        double Su = vonKarmanSu(X[2], freq);
        double sign = rx[i] / abs(rx[i]);
        // cout << rx[i] << " " << sign << endl;
        px[i] = sign * sqrt(1.0/nf*Su*df*(rx[i]*rx[i])/(1+rx[i]*rx[i]));
        qx[i] = sign * sqrt(1.0/nf*Su*df/(1+rx[i]*rx[i]));

        // Calculate py, qy
        double Sv = vonKarmanSv(X[2], freq);
        sign = ry[i] / abs(ry[i]);
        py[i] = sign * sqrt(1.0/nf*Sv*df*(ry[i]*ry[i])/(1+ry[i]*ry[i]));
        qy[i] = sign * sqrt(1.0/nf*Sv*df/(1+ry[i]*ry[i]));

        // Calculate pz, qz
        double Sw = vonKarmanSw(X[2], freq);
        sign = rz[i] / abs(rz[i]);
        pz[i] = sign * sqrt(1.0/nf*Sw*df*(rz[i]*rz[i])/(1+rz[i]*rz[i]));
        qz[i] = sign * sqrt(1.0/nf*Sw*df/(1+rz[i]*rz[i]));
        // cout << pz[i] << " " << qz[i] << endl;
    }

    vector<double> kx(nf);
    vector<double> ky(nf);
    vector<double> kz(nf);
    for(int i = 0; i < nf; ++i)
    {
        kx[i] = py[i]*qz[i] - pz[i]*qy[i];
        ky[i] = -(px[i]*qz[i] - pz[i]*qx[i]);
        kz[i] = px[i]*qy[i] - py[i]*qx[i];
        double mag = sqrt(kx[i]*kx[i]+ky[i]*ky[i]+kz[i]*kz[i]);
        kx[i] /= mag;
        ky[i] /= mag;
        kz[i] /= mag;

        // cout << kx[i] << " " << ky[i] << " " << kz[i] << endl;
    }

    vector<double> U(3, 0.0);

    double Uav = meanVelocity(X[2]);
    double gamma = tuningFactor(X[2]);
    double Cx = getCx();
    double Cy = getCy();
    double Cz = getCz();
    double Lsx = Uav / (gamma * Cx * freq);
    double Lsy = Uav / (gamma * Cy * freq);
    double Lsz = Uav / (gamma * Cz * freq);
    X[0] /= Lsx;
    X[1] /= Lsy;
    X[2] /= Lsz;

    for(int i = 0; i < nf; ++i)
    {
        double kx_ft = kx[i]*X[0] + ky[i]*X[1] + kz[i]*X[2] + 2*M_PI*fmn[i]*t;
        U[0] += px[i] * cos(kx_ft) + qx[i] * sin(kx_ft);
        U[1] += py[i] * cos(kx_ft) + qy[i] * sin(kx_ft);
        U[2] += pz[i] * cos(kx_ft) + qz[i] * sin(kx_ft);
    }

    cout << U[0] << " " << U[1] << " " << U[2] << endl;
    return U;
}
