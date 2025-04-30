#include "solver.h"
#include "phys_model.h"
#include "disk.h"
//#include "./lib/cosmos/astroconst.h"

using namespace phys_model::equations;

solver::TRungeKuttaSolver::TRungeKuttaSolver()
{

}

#ifdef STEADY_DRIFT
void solver::TRungeKuttaSolver::CalcCoefficients(const Tfloat_vector v0, const Tfloat dt, const Tfloat time)
{
    // radial distance from the star, au
    const double r_au = v0[0];
    // grain cross-section radius, cm
    const double a_cm = v0[1] * phys_model::scales::Get("a");
    // radial drift speed of the grain, non-dimensional (in units of scale v0)
    const double v_r = phys_model::equations::TerminalDriftSpeed(r_au, a_cm, disk::GetBeta_k());
    // non-dimensional azimuthal speed of grain
    double v_phi = disk::V_k(r_au) / phys_model::scales::Get("v");
    // specific angular momentum of the grain, non-dimensional (in units of scale j0)
    const double j = v_phi * r_au;

    /// Pointer to the solver coefficients
    Tfloat_vector_2d* cPtr = GetCoefficientsPtr();

    (*cPtr)[0][0] = Drdt(r_au, v_r, j, a_cm);
    (*cPtr)[0][1] = Dadt(r_au, v_r, j, a_cm);

    (*cPtr)[1][0] = Drdt(r_au + 0.5 * dt * GetC(0, 0), v_r, j, a_cm + 0.5 * dt * GetC(0, 1));
    (*cPtr)[1][1] = Dadt(r_au + 0.5 * dt * GetC(0, 0), v_r, j, a_cm + 0.5 * dt * GetC(0, 1));

    (*cPtr)[2][0] = Drdt(r_au + 0.5 * dt * GetC(1, 0), v_r, j, a_cm + 0.5 * dt * GetC(1, 1));
    (*cPtr)[2][1] = Dadt(r_au + 0.5 * dt * GetC(1, 0), v_r, j, a_cm + 0.5 * dt * GetC(1, 1));

    (*cPtr)[3][0] = Drdt(r_au + dt * GetC(2, 0), v_r, j, a_cm + dt * GetC(2, 1));
    (*cPtr)[3][1] = Dadt(r_au + dt * GetC(2, 0), v_r, j, a_cm + dt * GetC(2, 1));
}
#else
void solver::TRungeKuttaSolver::CalcCoefficients(const Tfloat_vector v0, const Tfloat dt, const Tfloat time)
{
    // radial distance from the star, au
    const double r_au = v0[0];
    // radial drift speed of the grain, non-dimensional (in units of scale v0)
    const double v_r = v0[1];
    // specific angular momentum of the grain, non-dimensional (in units of scale j0)
    const double j = v0[2];
    // grain cross-section radius, cm
    const double a_cm = v0[3] * phys_model::scales::Get("a");

    /// Pointer to the solver coefficients
    Tfloat_vector_2d* cPtr = GetCoefficientsPtr();

    (*cPtr)[0][0] = Drdt(r_au, v_r, j, a_cm);
    (*cPtr)[0][1] = Dvrdt(r_au, v_r, j, a_cm);
    (*cPtr)[0][2] = Djdt(r_au, v_r, j, a_cm);
    (*cPtr)[0][3] = Dadt(r_au, v_r, j, a_cm);

    (*cPtr)[1][0] = Drdt(r_au + 0.5 * dt * GetC(0, 0), v_r + 0.5 * dt * GetC(0, 1), j + 0.5 * dt * GetC(0, 2), a_cm + 0.5 * dt * GetC(0, 3));
    (*cPtr)[1][1] = Dvrdt(r_au + 0.5 * dt * GetC(0, 0), v_r + 0.5 * dt * GetC(0, 1), j + 0.5 * dt * GetC(0, 2), a_cm + 0.5 * dt * GetC(0, 3));
    (*cPtr)[1][2] = Djdt(r_au + 0.5 * dt * GetC(0, 0), v_r + 0.5 * dt * GetC(0, 1), j + 0.5 * dt * GetC(0, 2), a_cm + 0.5 * dt * GetC(0, 3));
    (*cPtr)[1][3] = Dadt(r_au + 0.5 * dt * GetC(0, 0), v_r + 0.5 * dt * GetC(0, 1), j + 0.5 * dt * GetC(0, 2), a_cm + 0.5 * dt * GetC(0, 3));

    (*cPtr)[2][0] = Drdt(r_au + 0.5 * dt * GetC(1, 0), v_r + 0.5 * dt * GetC(1, 1), j + 0.5 * dt * GetC(1, 2), a_cm + 0.5 * dt * GetC(1, 3));
    (*cPtr)[2][1] = Dvrdt(r_au + 0.5 * dt * GetC(1, 0), v_r + 0.5 * dt * GetC(1, 1), j + 0.5 * dt * GetC(1, 2), a_cm + 0.5 * dt * GetC(1, 3));
    (*cPtr)[2][2] = Djdt(r_au + 0.5 * dt * GetC(1, 0), v_r + 0.5 * dt * GetC(1, 1), j + 0.5 * dt * GetC(1, 2), a_cm + 0.5 * dt * GetC(1, 3));
    (*cPtr)[2][3] = Dadt(r_au + 0.5 * dt * GetC(1, 0), v_r + 0.5 * dt * GetC(1, 1), j + 0.5 * dt * GetC(1, 2), a_cm + 0.5 * dt * GetC(1, 3));


    (*cPtr)[3][0] = Drdt(r_au + dt * GetC(2, 0), v_r + dt * GetC(2, 1), j + dt * GetC(2, 2), a_cm + dt * GetC(2, 3));
    (*cPtr)[3][1] = Dvrdt(r_au + dt * GetC(2, 0), v_r + dt * GetC(2, 1), j + dt * GetC(2, 2), a_cm + dt * GetC(2, 3));
    (*cPtr)[3][2] = Djdt(r_au + dt * GetC(2, 0), v_r + dt * GetC(2, 1), j + dt * GetC(2, 2), a_cm + dt * GetC(2, 3));
    (*cPtr)[3][3] = Dadt(r_au + dt * GetC(2, 0), v_r + dt * GetC(2, 1), j + dt * GetC(2, 2), a_cm + dt * GetC(2, 3));
}
#endif
