#include "WindProfile.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    WindProfile wp;
    Info << wp.meanVelocity(0.364) << endl;
    Info << wp.turbIntensityU(0.364) << endl;
    Info << wp.turbIntensityV(0.364) << endl;
    Info << wp.turbIntensityW(0.364) << endl;
    Info << wp.turbLengthScaleU(0.364) << endl;
    Info << wp.turbLengthScaleV(0.364) << endl;
    Info << wp.turbLengthScaleW(0.364) << endl;
    Info << wp.vonKarmanSu(0.35, 1) << endl;
    Info << wp.gamma(0.364) << endl;

    return 0;
}
