#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED

#include "./lib/numeric/ode.h"

/// @brief Namespace - interface to the equations solver
namespace solver
{
    /// @brief Class - Runge-Kutta solver for the given equations system
    class TRungeKuttaSolver : public TBaseRungeKuttaSolver
	{
	private:
	public:
		/// @brief Constructor
		TRungeKuttaSolver();

		void CalcCoefficients(const Tfloat_vector v0, const Tfloat _dt, const Tfloat _time);
		void UserCheck() {};
	};
}

#endif // SOLVER_H_INCLUDED
