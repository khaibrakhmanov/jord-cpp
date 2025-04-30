#ifndef DISK_H_INCLUDED
#define DISK_H_INCLUDED

/// @brief Namespace - interface to the disk model
namespace disk
{
    /// @brief Make calculations to set up the disk model. Run before use the interface!
    void Prepare();

    /// @brief Set the degree of subkeplerian rotation from the input
    /// @param vphi2vk - degree defined as: vphi = vphi2vk * vk
    void SetSubkeplerDegree(const double vphi2vk);

    /// @brief Set disk model parameters from input
    /// @param alpha - turbulence parameter
    /// @param mdot - accretion rate, in units of 10^{-8} M_sun/year
    /// @param m - stellar mass, in units of solar masses
    /// @param vphi2vk - degree of subkeplerian rotation
    void SetParams(const double alpha, const double mdot, const double m, const double vphi2vk);

    /// @brief Set spatial disk boundaries in the r-direction
    /// @param r_in - inner radius of the disk
    /// @param r_out - outer radius of the disk
    void SetBoundaries(const double r_in, const double r_out);

    /// @brief Get degree of the subkeplerian rotation
    /// @return beta_k defined as beta_k = 2 * (1 - v_phi / v_k)
    double GetBeta_k();

    /// @brief Radial speed of the accretion flow
    /// @param r_au - radial distance from the star, au
    /// @return v_r [cm/s]
    double V_r(const double r_au);

    /// @brief Keplerian speed
    /// @param r_au - radial distance from the star, au
    /// @return v_k [cm/s]
    double V_k(const double r_au);

    /// @brief Azimuthal speed of the accretion flow
    /// @param r_au - radial distance from the star, au
    /// @return v_phi [cm/s]
    double V_phi(const double r_au);

    /// @brief Density of the matter in the disk
    /// @param r_au - radial distance from the star, au
    /// @return rho [g/cm^3]
    double Rho(const double r_au);

    /// @brief Thermal speed in the gas
    /// @param r_au - radial distance from the star, au
    /// @return v_th [cm/s]
    double V_th(const double r_au);
}

#endif // DISK_H_INCLUDED
