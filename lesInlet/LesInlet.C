#include "LesInlet.H"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>

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
    normal_distribution<double> dist_r(0.0,1.0);
    normal_distribution<double> dist_n(freq, df);

    auto gen_r = std::bind(dist_r, mt);
    vector<double> rx(nf);
    vector<double> ry(nf);
    vector<double> rz(nf);
    generate(begin(rx), end(rx), gen_r);
    generate(begin(ry), end(ry), gen_r);
    generate(begin(rz), end(rz), gen_r);

    // for (auto i : r)
    // {
    //     cout << i << " ";
    // }

    auto gen_fmn = std::bind(dist_n, mt);
    vector<double> fmn(nf);
    generate(begin(fmn), end(fmn), gen_fmn);

    vector<double> px(nf);
    vector<double> py(nf);
    vector<double> pz(nf);

    vector<double> qx(nf);
    vector<double> qy(nf);
    vector<double> qz(nf);

    for(int i = 0; i < nf; i++)
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
        cout << pz[i] << " " << qz[i] << endl;
    }

    vector<double> U(3);
    return U;
}
