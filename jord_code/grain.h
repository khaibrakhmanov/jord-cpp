#ifndef GRAIN_H_INCLUDED
#define GRAIN_H_INCLUDED

/// @brief Namespace - interface to the dust grain properties
namespace grain
{
    /// @brief Set spherical grain parameters from input
    /// @param a - cross-section radius, cm
    /// @param rho - density, g/cm^3
    /// @param Yd - dust-to-gas mass ratio
    void SetParameters(const double a, const double rho, const double Yd);

    /// @brief Get dust-to-gas mass ratio
    /// @return Yd
    double GetYd();

    /// @brief Get grain density
    /// @return rho, g/cm^3
    double GetRho();

    /// @brief Get grain cross-section radius
    /// @return a, cm
    double GetRadius();
}

#endif // GRAIN_H_INCLUDED
