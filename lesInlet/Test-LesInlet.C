#include "LesInlet.H"
#include "IOstreams.H"
#include <vector>

using namespace Foam;

int main()
{
    LesInlet wp;
    vector<double> X(3);
    X[0] = 0;
    X[1] = 0;
    X[2] = 0.364;
    double t = 0;
    vector<double> U = wp.Um(10, X, t);

    return 0;
}
