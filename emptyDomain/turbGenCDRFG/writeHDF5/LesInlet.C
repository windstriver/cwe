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

LesInlet::LesInlet(double deltaT)
{
    fmax = 1.0 / (5.0 * deltaT);
    nm = int(fmax / 2);
    nf = 100;
    fmin = fmax / nm / 2.0;
    df = (fmax - fmin) / (nm - 1);
}

double LesInlet::getFmax()
{
    return fmax;
}

std::vector<std::vector<std::vector<double> > >
LesInlet::UturbFreqSeg(double freq,
                       std::vector<std::vector<double> > xMat,
                       std::vector<double> tt)
{
    // Std::Vector to store coordinates of inlet point matrix
    std::vector<double> X(xMat.size(), 0.0);
    std::vector<double> Y(xMat.size(), 0.0);
    std::vector<double> Z(xMat.size(), 0.0);
    for (unsigned int i = 0; i < xMat.size(); ++i)
    {
        X[i] = xMat[i][0];
        Y[i] = xMat[i][1];
        Z[i] = xMat[i][2];
        // std::cout << X[i] << " " << Y[i] << " " << Z[i] << std::endl;
    }

    // Extract wind profile
    std::vector<double> Uav(Z.size(), 0.0);
    std::vector<double> gamma(Z.size(), 0.0);
    std::vector<double> Lsx(Z.size(), 0.0);
    std::vector<double> Lsy(Z.size(), 0.0);
    std::vector<double> Lsz(Z.size(), 0.0);
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
        // std::cout << X[nPt] << " " << Y[nPt] << " " << Z[nPt] << std::endl;
    }

    // Random generator
    std::random_device rd;    // non-deterministic generator
    std::mt19937 mt(rd());    // to seed mersenne twister
    // standard normal distribution
    std::normal_distribution<double> dist_r(0.0, 1.0);
    // normal distribution with mean of freq and std of df
    std::normal_distribution<double> dist_n(freq, df);

    // Std::Vector to store the standard normal r.v. ri
    std::vector<double> rx(nf);
    std::vector<double> ry(nf);
    std::vector<double> rz(nf);
    for(int i = 0; i < nf; ++i)
    {
        rx[i] = dist_r(mt);
        ry[i] = dist_r(mt);
        rz[i] = dist_r(mt);
    }

    // Std::Vector to store the normal distributed fmn ~ N(freq, df)
    std::vector<double> fmn(nf);
    for(int i = 0; i < nf; ++i)
    {
        fmn[i] = dist_n(mt);
    }

    // Std::Vector to calculate the p and q
    std::vector<std::vector<double> > px(Z.size(), std::vector<double>(nf));
    std::vector<std::vector<double> > py(Z.size(), std::vector<double>(nf));
    std::vector<std::vector<double> > pz(Z.size(), std::vector<double>(nf));
    std::vector<std::vector<double> > qx(Z.size(), std::vector<double>(nf));
    std::vector<std::vector<double> > qy(Z.size(), std::vector<double>(nf));
    std::vector<std::vector<double> > qz(Z.size(), std::vector<double>(nf));
    for (unsigned int nPt = 0; nPt < Z.size(); ++nPt)
    {
        // Calculate spectrum corresponding to point nPt
        double Su = vonKarmanSu(Z[nPt], freq);
        double Sv = vonKarmanSv(Z[nPt], freq);
        double Sw = vonKarmanSw(Z[nPt], freq);

        // std::cout << "von Karman spectrum:\n" << Su << " " << Sv << " " << Sw << std::endl;

        // Calculate p, q for different fmn around freq
        for(int i = 0; i < nf; ++i)
        {
            // Calculate px, qx
            double sign = rx[i] / std::abs(rx[i]);
            // std::cout << rx[i] << " " << sign << std::endl;
            px[nPt][i] = sign * std::sqrt(1.0/nf*Su*df*(rx[i]*rx[i])/(1+rx[i]*rx[i]));
            qx[nPt][i] = sign * std::sqrt(1.0/nf*Su*df/(1+rx[i]*rx[i]));

            // Calculate py, qy
            sign = ry[i] / std::abs(ry[i]);
            py[nPt][i] = sign *           std::sqrt(1.0/nf*Sv*df*(ry[i]*ry[i])/(1+ry[i]*ry[i]));
            qy[nPt][i] = sign * std::sqrt(1.0/nf*Sv*df/(1+ry[i]*ry[i]));

            // Calculate pz, qz
            sign = rz[i] / std::abs(rz[i]);
            pz[nPt][i] = sign * std::sqrt(1.0/nf*Sw*df*(rz[i]*rz[i])/(1+rz[i]*rz[i]));
            qz[nPt][i] = sign * std::sqrt(1.0/nf*Sw*df/(1+rz[i]*rz[i]));

            // std::cout << px[nPt][i] << " " << py[nPt][i] << " " << pz[nPt][i] << std::endl;
        }
    }

    // Calculate k to satisfy the divergence free condition
    std::vector<std::vector<double> > kx(Z.size(), std::vector<double>(nf));
    std::vector<std::vector<double> > ky(Z.size(), std::vector<double>(nf));
    std::vector<std::vector<double> > kz(Z.size(), std::vector<double>(nf));
    for (unsigned int nPt = 0; nPt < Z.size(); ++nPt)
    {
        for (int i = 0; i < nf; ++i)
        {
            kx[nPt][i] = py[nPt][i]*qz[nPt][i] - pz[nPt][i]*qy[nPt][i];
            ky[nPt][i] = -(px[nPt][i]*qz[nPt][i] - pz[nPt][i]*qx[nPt][i]);
            kz[nPt][i] = px[nPt][i]*qy[nPt][i] - py[nPt][i]*qx[nPt][i];
            double mag = std::sqrt(kx[nPt][i]*kx[nPt][i]+ky[nPt][i]*ky[nPt][i]+kz[nPt][i]*kz[nPt][i]);
            kx[nPt][i] /= mag;
            ky[nPt][i] /= mag;
            kz[nPt][i] /= mag;
        }
    }

    // Add for frequency samples around freq
    std::vector<std::vector<double> > Ux(Z.size(), std::vector<double>(tt.size()));
    std::vector<std::vector<double> > Uy(Z.size(), std::vector<double>(tt.size()));
    std::vector<std::vector<double> > Uz(Z.size(), std::vector<double>(tt.size()));
    for (unsigned nPt = 0; nPt < Z.size(); ++nPt)
    {
        for (unsigned ts = 0; ts < tt.size(); ++ts)
        {
            for(int i = 0; i < nf; ++i)
            {
                double kx_ft = kx[nPt][i]*X[nPt] + ky[nPt][i]*Y[nPt] + kz[nPt][i]*Z[nPt] + 2*M_PI*fmn[i]*tt[ts];
                Ux[nPt][ts] += px[nPt][i] * std::cos(kx_ft) + qx[nPt][i] * std::sin(kx_ft);
                Uy[nPt][ts] += py[nPt][i] * std::cos(kx_ft) + qy[nPt][i] * std::sin(kx_ft);
                Uz[nPt][ts] += pz[nPt][i] * std::cos(kx_ft) + qz[nPt][i] * std::sin(kx_ft);
            }
        }
    }

    std::vector<std::vector<std::vector<double> > >
        U(Z.size(),
          std::vector<std::vector<double> >(tt.size(), std::vector<double>(3)));

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

std::vector<std::vector<std::vector<double> > >
LesInlet::Uturb(std::vector<std::vector<double> > xMat,
                std::vector<double> tt)
{
    // Frequency std::vector
    std::vector<double> fm(nm, 0.0);
    for(int i = 0; i < nm; ++i)
    {
        fm[i] = fmin + i * df;
        // cout << fm[i] << endl;
    }

    // Add different frequency segments
    std::vector<std::vector<double> > Ux(xMat.size(), std::vector<double>(tt.size()));
    std::vector<std::vector<double> > Uy(xMat.size(), std::vector<double>(tt.size()));
    std::vector<std::vector<double> > Uz(xMat.size(), std::vector<double>(tt.size()));

    std::cout << "\nStart generating turbulence:" << std::endl;
    for(int i = 0; i < nm; ++i)
    {
        std::cout << "Frequency segment: " << fm[i] << "(Maximum freq.: " << fmax << ")" << std::endl;
        std::vector<std::vector<std::vector<double> > > Um= UturbFreqSeg(fm[i], xMat, tt);
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

    std::vector<std::vector<std::vector<double> > >
        U(xMat.size(),
          std::vector<std::vector<double> >(tt.size(), std::vector<double>(3)));

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

std::vector<double> LesInlet::Umean(std::vector<std::vector<double> > xMat)
{
    std::vector<double> Z(xMat.size(), 0.0);
    std::vector<double> Uav(Z.size(), 0.0);
    for (unsigned int i = 0; i < xMat.size(); i++)
    {
        Z[i] = xMat[i][2];
        Uav[i] = meanVelocity(Z[i]);
    }

    return Uav;
}
