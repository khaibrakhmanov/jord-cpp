#include "phys_model.h"
#include "grain.h"
#include "disk.h"
#include <cmath>
#include "./lib/coding/utypes.h"

namespace phys_model
{
    /// @brief Dictionary - model scales
    TFloatValueDict scalesDict{ {"r", 1.0}, {"v", 1.0}, {"j", 1.0}, {"t", 1.0}, {"a", 1.0}};

    namespace parameters
    {
        /// @brief The degree of subkeplerian rotation defined as: beta_k = 2 * (1 - v_phi / v_k)
        double beta_k = 0.01;
    }
}

void phys_model::scales::Set(const double r0, const double v0, const double a0)
{
    scalesDict["r"] = r0;
    scalesDict["v"] = v0;
    scalesDict["a"] = a0;

    scalesDict["t"] = r0 / v0;
    scalesDict["j"] = v0 * r0;
}

double phys_model::scales::Get(std::string name)
{
    return phys_model::scalesDict[name];
}


void phys_model::parameters::Set(const double beta_k)
{
    parameters::beta_k = beta_k;
};

double phys_model::equations::Dadt(const double r_au, const double vr, const double j, const double a_cm)
{
#ifdef COAGULATION
#ifdef STEADY_DRIFT
    double v_r = phys_model::equations::TerminalDriftSpeed(r_au, a_cm, disk::GetBeta_k());
#endif
    // non-dimensional radial speed of gas
    double u_r = disk::V_r(r_au) / scales::Get("v");
    // non-dimensional azimuthal speed of gas
    double u_phi = disk::V_phi(r_au) / scales::Get("v");
    // non-dimensional azimuthal speed of grain
    double v_phi = j / r_au;
#ifdef STEADY_DRIFT
    // total (relative to gas) speed of the grain in the (r-z)-plane
    double v = sqrt(pow(v_r - u_r, 2) + pow(v_phi - u_phi, 2));
#else
    // total (relative to gas) speed of the grain in the (r-z)-plane
    double v = sqrt(pow(vr - u_r, 2) + pow(v_phi - u_phi, 2));
#endif

    return scales::Get("r") * grain::GetYd() * disk::Rho(r_au) / (scales::Get("a") * 4 * grain::GetRho()) * v;
#else
    return 0.0;
#endif
}

double phys_model::equations::TerminalDriftSpeed(const double r_au, const double a_cm, const double beta_k)
{
    // non-dimensional radial speed of gas
    double u_r = disk::V_r(r_au) / scales::Get("v");
    // non-dimensional stopping time at current state
    double tau_s = coefficients::StoppingTime(r_au, a_cm) / scales::Get("t");

    return (u_r - beta_k * tau_s / pow(r_au, 2.0)) / (1.0 + pow(tau_s, 2.0) / pow(r_au, 3.0));
}

double phys_model::equations::Drdt(const double r_au, const double vr, const double j, const double a_cm)
{
#ifdef STEADY_DRIFT
    return TerminalDriftSpeed(r_au, a_cm, disk::GetBeta_k());
#else
    return vr;
#endif
}

#ifndef STEADY_DRIFT
double phys_model::equations::Dvrdt(const double r_au, const double vr, const double j, const double a_cm)
{
    // cube of the radial coordinate
    double rrr = pow(r_au, 3.0);
    // non-dimensional stopping time at current state
    double tau_s = coefficients::StoppingTime(r_au, a_cm) / scales::Get("t");
    // non-dimensional radial speed of gas
    double u_r = disk::V_r(r_au) / scales::Get("v");

    return pow(j, 2.0) / rrr - pow(r_au, -2.0) - (vr - u_r) / tau_s;
}

double phys_model::equations::Djdt(const double r_au, const double vr, const double j, const double a_cm)
{
    // non-dimensional stopping time at current state
    double tau_s = coefficients::StoppingTime(r_au, a_cm) / scales::Get("t");
    // non-dimensional azimuthal speed of grain
    double v_phi = j / r_au;
    // non-dimensional azimuthal speed of gas
    double u_phi = disk::V_phi(r_au) / scales::Get("v");

    return (-r_au / tau_s) * (v_phi - u_phi);
}
#endif

double phys_model::coefficients::StoppingTime(const double r_au, const double a_cm)
{
         // the solution for the case of Epstein drag
    if ((a_cm / scales::Get("a")) < 1.5 * disk::lambda(r_au))
    {
        return grain::GetRho() * a_cm / (disk::Rho(r_au) * disk::V_th(r_au)); 
    }
    else // the solution for the case of Stokes drag
    {
        return 2 * grain::GetRho() * pow(a_cm, 2) / (3 * disk::Rho(r_au) * disk::V_th(r_au) * disk::lambda(r_au)); 
    }
}
