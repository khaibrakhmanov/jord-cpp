#include "disk.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include "./lib/cosmos/astroconst.h"
#include "./lib/physics/physconst.h"

using namespace astro_const;
using namespace physical_const;

namespace disk
{
    /// gas molecular weight
    double mu = 2.3;

    /// turbulence parameter
    double alpha = 0.01;
    /// accretion rate, in units of 10^{-8} M_sun/year
    double m_dot = 1.0;
    /// stellar mass, in units of solar masses
    double m = 1.0;

    /// @brief degree of subkeplerian rotation defined as: vphi = vphi2vk * vk
    double vphi2vk = 1.0;
    /// @brief degree of the subkeplerian rotation defined as: beta_k = 2 * (1 - v_phi / v_k)
    double beta_k = 0.01;

    /// @brief inner radius of the disk, au
    double r_in = 0.05;
    /// @brief outer radius of the disk, au
    double r_out = 100.0;

    // ----- constants in the solution of the model equations -----

    double C_ur = 1.0;
    double C_vk = 1.0;
    double C_rho = 1.0;
    double C_vth = 1.0;
    double C_T = 1.0;
    // ------------------------------------------------------------
}

void disk::SetParams(const double alpha, const double mdot, const double m, const double vphi2vk)
{
    disk::alpha = alpha;
    disk::m_dot = mdot;
    disk::m = m;
    disk::vphi2vk = vphi2vk;
    disk::beta_k = 2.0 * (1.0 - disk::vphi2vk);
}

void disk::SetBoundaries(const double in, const double out)
{
    r_in = in;
    r_out = out;
};

double disk::GetBeta_k()
{
    return beta_k;
}
void disk::Prepare()
{
    C_ur = 30 * pow(alpha / 0.01, 0.75) * sqrt(m_dot) * pow(m, -0.125);
    C_vk = sqrt(G * M_sun /au);
    C_rho = 2.5e-10 * pow(disk::m, 0.4375) * pow(disk::alpha / 0.01, -0.625) * pow(disk::m_dot, 0.25);
    C_vth = sqrt(8 * k / (M_PI * mu * m_p));
    C_T = 240 * pow(disk::m, 0.375) * pow(disk::alpha / 0.01, -0.25) * sqrt(disk::m_dot);
}

void disk::SetSubkeplerDegree(const double arg)
{
    vphi2vk = arg;
}

double disk::V_r(const double r_au)
{
    return -C_ur * pow(r_au, -0.625);
}

double disk::V_k(const double r_au)
{
    return C_vk * sqrt(disk::m / r_au);
}

double disk::V_phi(const double r_au)
{
    return disk::vphi2vk * V_k(r_au);
}

double disk::Rho(const double r_au)
{
    return C_rho * pow(r_au, -1.3125);
}

double T(const double r_au)
{
    return disk::C_T * pow(r_au, -1.125);
}

double disk::V_th(const double r_au)
{
    return C_vth * sqrt(T(r_au));
}
