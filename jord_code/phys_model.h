#ifndef EQS_H_INCLUDED
#define EQS_H_INCLUDED

#include <string>

//#define COAGULATION
//#define STEADY_DRIFT

/// @brief Namespace - interface to the physical model (equations, coefficients, scales, etc.)
namespace phys_model
{
    /// @brief Namespace - interface to model scales
    namespace scales
    {
        /// @brief Set the basic scales from input
        /// @param r0 - coordinate scale, cm
        /// @param v0 - veloctity scale, cm/s
        /// @param a0 - cross-section radius scale, cm
        void Set(const double r0, const double v0, const double a0);

        /// @brief Get the scale of a given value
        /// @param name - value name
        /// @return scale, CGS
        double Get(std::string name);
    }

    /// @brief Namespace - interface to model parameters
    namespace parameters
    {
        /// @brief Set the degree of subkeplerian rotation
        /// @param beta_k - the degree of subkeplerian rotation defined as: beta_k = 2 * (1 - v_phi / v_k)
        void Set(const double beta_k);
    }

    /// @brief Namespace - interface to the coefficients of model equations
    namespace coefficients
    {
        /// @brief Stopping time
        /// @param r_au - radial distance from the star, au
        /// @param a_cm - grain cross-section radius, cm
        /// @return t_stop, s
        double StoppingTime(const double r_au, const double a_cm);
    }

    /// @brief Namespace - interface to the model equations
    namespace equations
    {
        /// @brief Derivative dr/dt, i.e. the right-hand-side of the equation for the r-coordinate
        /// @param r_au - radial distance from the star, au
        /// @param vr - radial drift speed of the grain, non-dimensional (in units of scale v0)
        /// @param j - specific angular momentum of the grain, non-dimensional (in units of scale j0)
        /// @param a_cm - grain cross-section radius, cm
        /// @return dr/dt, non-dimensional
        double Drdt(const double r_au, const double vr, const double j, const double a_cm);

#ifndef STEADY_DRIFT
        /// @brief Derivative d(v_r)/dt, i.e. the right-hand-side of the equation for the radial drift speed v_r
        /// @param r_au - radial distance from the star, au
        /// @param vr - radial drift speed of the grain, non-dimensional (in units of scale v0)
        /// @param j - specific angular momentum of the grain, non-dimensional (in units of scale j0)
        /// @param a_cm - grain cross-section radius, cm
        /// @return d(v_r)/dt, non-dimensional
        double Dvrdt(const double r_au, const double vr, const double j, const double a_cm);

        /// @brief Derivative dj/dt, i.e. the right-hand-side of the equation for the specific angular momentum
        /// @param r_au - radial distance from the star, au
        /// @param vr - radial drift speed of the grain, non-dimensional (in units of scale v0)
        /// @param j - specific angular momentum of the grain, non-dimensional (in units of scale j0)
        /// @param a_cm - grain cross-section radius, cm
        /// @return dj/dt, non-dimensional
        double Djdt(const double r_au, const double vr, const double j, const double a_cm);
#endif
        /// @brief Derivative da/dt, i.e. the right-hand-side of the equation for the grain cross-section radius
        /// @param r_au - radial distance from the star, au
        /// @param vr - radial drift speed of the grain, non-dimensional (in units of scale v0)
        /// @param j - specific angular momentum of the grain, non-dimensional (in units of scale j0)
        /// @param a_cm - grain cross-section radius, cm
        /// @return da/dt, non-dimensional
        double Dadt(const double r_au, const double vr, const double j, const double a_cm);

        /// @brief Terminal radial drift speed
        /// @param r_au - radial distance from the star, au
        /// @param a_cm - grain cross-section radius, cm
        /// @param beta_k - the degree of subkeplerian rotation defined as: beta_k = 2 * (1 - v_phi / v_k)
        /// @return v_r, non-dimensional (in units of scale v0)
        double TerminalDriftSpeed(const double r_au, const double a_cm, const double beta_k);
    }
}

#endif // EQS_H_INCLUDED
