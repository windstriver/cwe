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

double LesInlet::getFmax()
{
    return fmax;
}

vector<vector<vector<double> > >
LesInlet::UturbFreqSeg(double freq,
                       vector<vector<double> > xMat,
                       vector<double> tt)
{
    // Vector to store coordinates of inlet point matrix
    vector<double> X(xMat.size(), 0.0);
    vector<double> Y(xMat.size(), 0.0);
    vector<double> Z(xMat.size(), 0.0);
    for (unsigned int i = 0; i < xMat.size(); ++i)
    {
        X[i] = xMat[i][0];
        Y[i] = xMat[i][1];
        Z[i] = xMat[i][2];
        // cout << X[i] << " " << Y[i] << " " << Z[i] << endl;
    }

    // Extract wind profile
    vector<double> Uav(Z.size(), 0.0);
    vector<double> gamma(Z.size(), 0.0);
    vector<double> Lsx(Z.size(), 0.0);
    vector<double> Lsy(Z.size(), 0.0);
    vector<double> Lsz(Z.size(), 0.0);
    double Cx = getCx();
    double Cy = getCy();
    double Cz = getCz();
    for (unsigned int nPt = 0; nPt < Z.size(); ++nPt)
    {
        Uav[nPt] = meanVelocity(Z[nPt]);
        // cout << Uav[nPt] << endl;
        gamma[nPt] = tuningFactor(Z[nPt]);
        // cout << gamma[nPt] << endl;
        Lsx[nPt] = Uav[nPt] / (gamma[nPt] * Cx * freq);
        Lsy[nPt] = Uav[nPt] / (gamma[nPt] * Cy * freq);
        Lsz[nPt] = Uav[nPt] / (gamma[nPt] * Cz * freq);
        // cout << Lsx[nPt] << " " << Lsy[nPt] << " " << Lsz[nPt] << endl;
        // scale the inlet points c.s.
        X[nPt] /= Lsx[nPt];
        Y[nPt] /= Lsy[nPt];
        Z[nPt] /= Lsz[nPt];
        // cout << X[nPt] << " " << Y[nPt] << " " << Z[nPt] << endl;
    }

    // Random generator
    random_device rd;    // non-deterministic generator
    mt19937 mt(rd());    // to seed mersenne twister
    // standard normal distribution
    normal_distribution<double> dist_r(0.0, 1.0);
    // normal distribution with mean of freq and std of df
    normal_distribution<double> dist_n(freq, df);

    // Vector to store the standard normal r.v. ri
    vector<double> rx(nf);
    vector<double> ry(nf);
    vector<double> rz(nf);
    for(int i = 0; i < nf; ++i)
    {
        rx[i] = dist_r(mt);
        ry[i] = dist_r(mt);
        rz[i] = dist_r(mt);
    }

    // Vector to store the normal distributed fmn ~ N(freq, df)
    vector<double> fmn(nf);
    for(int i = 0; i < nf; ++i)
    {
        fmn[i] = dist_n(mt);
    }

    // Vector to calculate the p and q
    vector<vector<double> > px(Z.size(), vector<double>(nf));
    vector<vector<double> > py(Z.size(), vector<double>(nf));
    vector<vector<double> > pz(Z.size(), vector<double>(nf));
    vector<vector<double> > qx(Z.size(), vector<double>(nf));
    vector<vector<double> > qy(Z.size(), vector<double>(nf));
    vector<vector<double> > qz(Z.size(), vector<double>(nf));
    for (unsigned int nPt = 0; nPt < Z.size(); ++nPt)
    {
        // Calculate spectrum corresponding to point nPt
        double Su = vonKarmanSu(Z[nPt], freq);
        double Sv = vonKarmanSv(Z[nPt], freq);
        double Sw = vonKarmanSw(Z[nPt], freq);

        // Calculate p, q for different fmn around freq
        for(int i = 0; i < nf; ++i)
        {
            // Calculate px, qx
            double sign = rx[i] / abs(rx[i]);
            // cout << rx[i] << " " << sign << endl;
            px[nPt][i] = sign * sqrt(1.0/nf*Su*df*(rx[i]*rx[i])/(1+rx[i]*rx[i]));
            qx[nPt][i] = sign * sqrt(1.0/nf*Su*df/(1+rx[i]*rx[i]));

            // Calculate py, qy
            sign = ry[i] / abs(ry[i]);
            py[nPt][i] = sign *           sqrt(1.0/nf*Sv*df*(ry[i]*ry[i])/(1+ry[i]*ry[i]));
            qy[nPt][i] = sign * sqrt(1.0/nf*Sv*df/(1+ry[i]*ry[i]));

            // Calculate pz, qz
            sign = rz[i] / abs(rz[i]);
            pz[nPt][i] = sign * sqrt(1.0/nf*Sw*df*(rz[i]*rz[i])/(1+rz[i]*rz[i]));
            qz[nPt][i] = sign * sqrt(1.0/nf*Sw*df/(1+rz[i]*rz[i]));
        }
    }

    // Calculate k to satisfy the divergence free condition
    vector<vector<double> > kx(Z.size(), vector<double>(nf));
    vector<vector<double> > ky(Z.size(), vector<double>(nf));
    vector<vector<double> > kz(Z.size(), vector<double>(nf));
    for (unsigned int nPt = 0; nPt < Z.size(); ++nPt)
    {
        for (int i = 0; i < nf; ++i)
        {
            kx[nPt][i] = py[nPt][i]*qz[nPt][i] - pz[nPt][i]*qy[nPt][i];
            ky[nPt][i] = -(px[nPt][i]*qz[nPt][i] - pz[nPt][i]*qx[nPt][i]);
            kz[nPt][i] = px[nPt][i]*qy[nPt][i] - py[nPt][i]*qx[nPt][i];
            double mag = sqrt(kx[nPt][i]*kx[nPt][i]+ky[nPt][i]*ky[nPt][i]+kz[nPt][i]*kz[nPt][i]);
            kx[nPt][i] /= mag;
            ky[nPt][i] /= mag;
            kz[nPt][i] /= mag;
        }
    }

    // Add for frequency samples around freq
    vector<vector<double> > Ux(Z.size(), vector<double>(tt.size()));
    vector<vector<double> > Uy(Z.size(), vector<double>(tt.size()));
    vector<vector<double> > Uz(Z.size(), vector<double>(tt.size()));
    for (unsigned nPt = 0; nPt < Z.size(); ++nPt)
    {
        for (unsigned ts = 0; ts < tt.size(); ++ts)
        {
            for(int i = 0; i < nf; ++i)
            {
                double kx_ft = kx[nPt][i]*X[nPt] + ky[nPt][i]*Y[nPt] + kz[nPt][i]*Z[nPt] + 2*M_PI*fmn[i]*tt[ts];
                Ux[nPt][ts] += px[nPt][i] * cos(kx_ft) + qx[nPt][i] * sin(kx_ft);
                Uy[nPt][ts] += py[nPt][i] * cos(kx_ft) + qy[nPt][i] * sin(kx_ft);
                Uz[nPt][ts] += pz[nPt][i] * cos(kx_ft) + qz[nPt][i] * sin(kx_ft);
            }
        }
    }

    vector<vector<vector<double> > >
        U(Z.size(),
          vector<vector<double> >(tt.size(), vector<double>(3)));

    for (unsigned int nPt = 0; nPt < Z.size(); ++nPt)
    {
        for (unsigned int ts = 0; ts < tt.size(); ++ts)
        {
            U[nPt][ts][0] = Ux[nPt][ts];
            U[nPt][ts][1] = Uy[nPt][ts];
            U[nPt][ts][2] = Uz[nPt][ts];
        }
    }

    return U;
}

vector<vector<vector<double> > >
LesInlet::Uturb(vector<vector<double> > xMat,
                vector<double> tt)
{
    // Frequency vector
    vector<double> fm(nm, 0.0);
    for(int i = 0; i < nm; ++i)
    {
        fm[i] = fmin + i * df;
        // cout << fm[i] << endl;
    }

    // Add different frequency segments
    vector<vector<double> > Ux(xMat.size(), vector<double>(tt.size()));
    vector<vector<double> > Uy(xMat.size(), vector<double>(tt.size()));
    vector<vector<double> > Uz(xMat.size(), vector<double>(tt.size()));

    for(int i = 0; i < nm; ++i)
    {
        vector<vector<vector<double> > > Um= UturbFreqSeg(fm[i], xMat, tt);
        for (unsigned nPt = 0; nPt < xMat.size(); ++nPt)
        {
            for (unsigned ts = 0; ts < tt.size(); ++ts)
            {
                Ux[nPt][ts] += Um[nPt][ts][0];
                Uy[nPt][ts] += Um[nPt][ts][1];
                Uz[nPt][ts] += Um[nPt][ts][2];
            }
        }
    }

    vector<vector<vector<double> > >
        U(xMat.size(),
          vector<vector<double> >(tt.size(), vector<double>(3)));

    for (unsigned int nPt = 0; nPt < xMat.size(); ++nPt)
    {
        for (unsigned int ts = 0; ts < tt.size(); ++ts)
        {
            U[nPt][ts][0] = Ux[nPt][ts];
            U[nPt][ts][1] = Uy[nPt][ts];
            U[nPt][ts][2] = Uz[nPt][ts];
        }
    }

    return U;
}

