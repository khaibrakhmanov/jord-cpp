#include "grain.h"
#include "disk.h"
#include "./lib/cosmos/astroconst.h"
#include <cmath>

using namespace astro_const;

namespace grain
{
    /// @brief grain density, g/cm^3
    double rho = 3.0;
    /// @brief grain cross-section radius, cm
    double a = 1e-5;
    /// @brief dust-to-gas mass ratio
    double massFraction = 0.01;
}

void grain::SetParameters(const double a, const double rho, const double Yd)
{
    grain::rho = rho;
    grain::a = a;
    grain::massFraction = Yd;
};

double grain::GetYd()
{
    return massFraction;
}

double grain::GetRho()
{
    return rho;
}

double grain::GetRadius()
{
    return a;
}