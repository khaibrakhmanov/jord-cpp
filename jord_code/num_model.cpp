#include "num_model.h"
#include <cmath>
#include "grain.h"
#include "disk.h"
#include "./lib/coding/params_reader.h"
#include <iomanip>
#include <iostream>
#include "solver.h"
#include <fstream>
#include "phys_model.h"
#include "./lib/cosmos/astroconst.h"
#include <algorithm>
#include <chrono>

namespace model
{
    /// Namespace - model parameters
    namespace parameters
    {
        /// current time, non-dimensional
        double t = 0.0;
        /// index of current time step
        std::size_t step = 0;
        /// time to stop simulation, non-dimensional
        double stopTime = 1.0;
        /// duration of the run, sec
        std::chrono::duration<double> runTime;
        /// current time step, non-dimensional
        double dt = 0.1;
        /// period for saving data
        std::size_t savePeriod = 1;
        /// period for screen output
        std::size_t outputPeriod = 1;
#ifdef STEADY_DRIFT
        /// number of variables
        const std::size_t varN = 2;
#else
        /// number of variables
        const std::size_t varN = 4;
#endif
        /// Print model parameters on screen
        void Print();
    }

    /// Namespace - model variables
    namespace variables
    {
        /// current radial coordinate of the grain, au
        double r_au = 1.0;
        /// current radial speed of the grain, non-dimensional
        double v_r = 1.0;
        /// current specific angular momentum of the grain, non-dimensional
        double j = 1.0;
        /// current cross-section radius of the grain, non-dimensional
        double a_nd = 1.0;
        /// vector of the solution of model equations at current time step
        Tfloat_vector sol;
    }

    /// Namespace - interface to model input/output functions
    namespace io
    {
        /// @brief Reader of model parameters
        TParameterReader reader;
        /// @brief Name of the file to save the simulation results
        std::ofstream outfile;
        /// @brief Open output data file
        void OpenOutputDataFile();
        /// @brief Save the simulation results at current time step
        void SaveCurrentData();
        /// @brief Print model summary to the summary file
        void PrintSummary();
    }

    /// Namespace - interface to the numerical scheme
    namespace scheme
    {
        /// @brief Solver to solve the model equations
        solver::TRungeKuttaSolver solver;
        /// @brief Set solver parameters
        void ConfigureSolver();
    }

    /// Namespace - interface to the steady state solution of the model equations
    namespace steady_solution
    {
        /// array of the radial grid coordinates
        Tfloat_vector rs;
        /// non-dimensional stopping time at the starting point
        double tauStop0;
        /// minimum non-dimensional stopping time on the radial grid
        double tauStopMin;
        /// array of non-dimensional stopping time values on the radial grid
        Tfloat_vector tauStops;
    }

    /// @brief Make analytical estimates to configure the scheme
    void AnalyticEstimates();
    /// @brief Prepare model: configure all units and set all parameters
    void Prepare();
    /// @brief Finalize run
    void Finalize();

}

void model::ReadParameters(std::string iniFile)
{
    try
    {
        io::reader.ReadParameters(iniFile);
    } catch (const std::runtime_error& e)
    {
        std::cerr << "Error while reading params file! Program will not run simulation correctly " << e.what() << std::endl;
    }
}

void model::scheme::ConfigureSolver()
{
    using namespace io;

    // We use Runge-Kutta schmeme of the 4th order of accuracy
    // with step size control

    solver.SetTime(0.0);
    solver.SetTimeStep(reader.GetDouble("dt0"));
    solver.SetOrder(4);
    solver.SetTolerance(reader.GetDouble("tol"));
    solver.SetVariablesNumber(parameters::varN);
    solver.ActivateStepSizeControl();

    // Min and max time steps are defined in comparison with
    // the minimum value of the stopping time in the disk, tauStopMin.
    // Time step should not be larger than the tauStopMin,
    // otherwise the graing trajectory will noy be modeled with sufficient accuracy

    double minStep = reader.GetDouble("dt_min") * steady_solution::tauStopMin;
    solver.SetMinStep(minStep);
    double maxStep = reader.GetDouble("dt_max") * steady_solution::tauStopMin;
    solver.SetMaxStep(maxStep);

    solver.AllocateSolution();
    solver.AllocateCoefficients();
}

void model::Finalize()
{
    // close output data file
    io::outfile.close();
    // print model summary to file
    io::PrintSummary();
}

void model::Prepare()
{
    using namespace io;

    // ----- set disk parameters -----

    disk::SetParams(reader.GetDouble("alpha"), reader.GetDouble("mdot"), reader.GetDouble("m"), reader.GetDouble("vphi2vk"));
    disk::SetBoundaries(reader.GetDouble("r_in"), reader.GetDouble("r_out"));
    disk::Prepare();

    // ----- set grain parameters -----

    grain::SetParameters(reader.GetDouble("a_0"), reader.GetDouble("rho"), reader.GetDouble("Yd"));

    // ----- Set physical model and equation parameters ------------
    //       Use properties of the Keplerian orbit at r = 1 au
    //       as scales of the coordinate and speed.
    //       The scale of the grain cross-section radius - its initial value
    phys_model::scales::Set(astro_const::au, disk::V_k(1.0), grain::GetRadius());

    // ----- set numerical model  parameters -----

    parameters::stopTime = reader.GetDouble("t_stop");
    parameters::savePeriod = reader.GetInt("save_period");
    parameters::outputPeriod = reader.GetInt("output_period");
    AllocateVector1D(&variables::sol, parameters::varN);
    AllocateVector1D(&steady_solution::rs, 100);
    AllocateVector1D(&steady_solution::tauStops, 100);
    AnalyticEstimates();

    // ----- set scheme parameters -----

    scheme::ConfigureSolver();

    // ----- prepare input/output interface -----

    try {
        OpenOutputDataFile();
    }
    catch (const std::runtime_error& e) {
        std::cerr << "Error! " << e.what() << std::endl;
    }
    outfile << "# t" << " " << "r" << " " << "v_r" << " " << "j" << " " << "a" << " " << "tau_s" << " " << "v_steady" <<  std::endl;
    outfile.precision(3);
    outfile << std::scientific;
}

void model::AnalyticEstimates()
{
    using namespace steady_solution;

    // stopping time at the initial point, s
    steady_solution::tauStop0 = phys_model::coefficients::StoppingTime(io::reader.GetDouble("r_0"), io::reader.GetDouble("a_0"));
    // stopping time at the initial point, non-dimensional
    steady_solution::tauStop0 = steady_solution::tauStop0 / phys_model::scales::Get("t");

    // radial grid cell size
    double dr = (io::reader.GetDouble("r_out") - io::reader.GetDouble("r_in")) / (100-1);
    // initial point in the radial grid
    steady_solution::rs[0] = io::reader.GetDouble("r_in");
    // stopping time at r[0]
    steady_solution::tauStops[0] = phys_model::coefficients::StoppingTime(rs[0], io::reader.GetDouble("a_0"));
    // calculate stopping time at each grid cell
    for (size_t i = 1; i < 100; i++)
    {
        steady_solution::rs[i] = steady_solution::rs[i - 1] + dr;
        steady_solution::tauStops[i] = phys_model::coefficients::StoppingTime(rs[i], io::reader.GetDouble("a_0"));
    }

    Tfloat_vector::iterator minValue = std::min_element(tauStops.begin(), tauStops.end());
    // minimum value of the stopping time on the grid, s
    steady_solution::tauStopMin = *minValue;
    // minimum value of the stopping time on the grid, non-dimensional
    steady_solution::tauStopMin = steady_solution::tauStopMin / phys_model::scales::Get("t");
}

void model::parameters::Print()
{
    std::cout << " Disk parameters: " << std::endl;
    std::cout << "  beta_k = " << disk::GetBeta_k() << std::endl;

    std::cout << " Scales: " << std::endl;
    std::cout << std::scientific;
    std::cout << std::setprecision(2);
    std::cout << "  t0       = " << phys_model::scales::Get("t") << " [s] = " << phys_model::scales::Get("t") / astro_const::year << " [yr]" << std::endl;
    std::cout << "  r0       = " << phys_model::scales::Get("r") / astro_const::au << " [au]" << std::endl;
    std::cout << "  v0       = " << phys_model::scales::Get("v") / 1e5 << " [km/s]" << std::endl;
    std::cout << "  j0       = " << phys_model::scales::Get("j") << " [cm^2/s]" << std::endl;
    std::cout << "  a0       = " << phys_model::scales::Get("a") << " [cm]" << std::endl;
    std::cout << "  tauStop (r0)  = " << steady_solution::tauStop0 << " [t0] = " << steady_solution::tauStop0 * phys_model::scales::Get("t") << " [s] = " << steady_solution::tauStop0 * phys_model::scales::Get("t") / astro_const::year << " [yr]" << std::endl;
    std::cout << "  tauStop (min) = " << steady_solution::tauStopMin << " [t0] = " << steady_solution::tauStopMin * phys_model::scales::Get("t") << " [s] = " << steady_solution::tauStopMin * phys_model::scales::Get("t") / astro_const::year << " [yr]" << std::endl;

    std::cout << " Initial state: " << std::endl;
    std::cout << "  r_0    = " << variables::r_au << " [au]" << std::endl;
    std::cout << "  v_0    = " << variables::v_r << " [v0]" << std::endl;
    std::cout << "  j_0    = " << variables::j << " [j0]" << std::endl;
    std::cout << "  a_0    = " << variables::a_nd * phys_model::scales::Get("a") << " [cm]" << std::endl;

    std::cout << " Grain parameters: " << std::endl;
    std::cout << "  rho    = " << grain::GetRho() << " [g/cm^3]" << std::endl;
    std::cout << "  Y_d    = " << grain::GetYd() << std::endl;

    std::cout << " Scheme parameters: " << std::endl;
    std::cout << "  dt0    = " << scheme::solver.GetTimeStep() << " [t0]" << std::endl;
    std::cout << "  dt_min = " << scheme::solver.GetMinStep() << " [t0]" << std::endl;
    std::cout << "  dt_max = " << scheme::solver.GetMaxStep() << " [t0]" << std::endl;

    std::cout << std::fixed;
}

void model::io::OpenOutputDataFile()
{
    outfile.open("data.dat", std::ofstream::out);
    if (!outfile.is_open()) {
        throw std::runtime_error("Could not open output data file!");
    }
}

void model::io::SaveCurrentData()
{
    using namespace variables;

    // current grain cross-section radius, cm
    double a_cm = a_nd * phys_model::scales::Get("a");

    // data are saved in the following order: r, v_r (radial drift speed), j, a, tau_stop, terminal v_r (steady-state solution)
    // all values are non-dimensional

    outfile << parameters::t << " " << r_au << " " << v_r << " " << j << " " << a_nd << " " << phys_model::coefficients::StoppingTime(sol[0], a_cm) << " " << phys_model::equations::TerminalDriftSpeed(r_au, a_cm, disk::GetBeta_k()) << std::endl;
}

void model::io::PrintSummary()
{
    std::ofstream summaryFile("model_summary.txt");
    if (summaryFile.is_open()) {
        summaryFile << "----- Scales -----" << std::endl;
        summaryFile << std::scientific;
        summaryFile << std::setprecision(2);
        summaryFile << "  t0       = " << phys_model::scales::Get("t") << " [s] = " << phys_model::scales::Get("t") / astro_const::year << " [yr]" << std::endl;
        summaryFile << "  r0       = " << phys_model::scales::Get("r") / astro_const::au << " [au]" << std::endl;
        summaryFile << "  v0       = " << phys_model::scales::Get("v") / 1e5 << " [km/s]" << std::endl;
        summaryFile << "  j0       = " << phys_model::scales::Get("j") << " [cm^2/s]" << std::endl;
        summaryFile << "  a0       = " << phys_model::scales::Get("a") << " [cm]" << std::endl;
        summaryFile << "  tauStop (r0)  = " << steady_solution::tauStop0 << " [t0] = " << steady_solution::tauStop0 * phys_model::scales::Get("t") << " [s] = " << steady_solution::tauStop0 * phys_model::scales::Get("t") / astro_const::year << " [yr]" << std::endl;
        summaryFile << "  tauStop (min) = " << steady_solution::tauStopMin << " [t0] = " << steady_solution::tauStopMin * phys_model::scales::Get("t") << " [s] = " << steady_solution::tauStopMin * phys_model::scales::Get("t") / astro_const::year << " [yr]" << std::endl;

        summaryFile << std::endl << "----- Initial state -----" << std::endl;
        summaryFile << "  r_0    = " << reader.GetDouble("r_0") << " [au]" << std::endl;
        summaryFile << "  v_0    = " << reader.GetDouble("v_0") << " [v0]" << std::endl;
        summaryFile << "  j_0    = " << reader.GetDouble("j_0") << " [j0]" << std::endl;
        summaryFile << "  a_0    = " << reader.GetDouble("a_0") << " [cm]" << std::endl;

        summaryFile << std::endl << "----- Grain parameters -----" << std::endl;
        summaryFile << "  rho    = " << grain::GetRho() << " [g/cm^3]" << std::endl;
        summaryFile << "  Y_d    = " << grain::GetYd() << std::endl;

        summaryFile << std::endl << "----- Scheme parameters -----" << std::endl;
        summaryFile << "  dt0    = " << scheme::solver.GetTimeStep() << " [t0]" << std::endl;
        summaryFile << "  dt_min = " << scheme::solver.GetMinStep() << " [t0]" << std::endl;
        summaryFile << "  dt_max = " << scheme::solver.GetMaxStep() << " [t0]" << std::endl;

        summaryFile << std::endl << "----- Run parameters -----" << std::endl;
        summaryFile << "  t_end    = " << parameters::t << " [t0] = "
            << parameters::t * phys_model::scales::Get("t") / astro_const::year << " [yr] " << std::endl;
        summaryFile << "  steps    = " << parameters::step << std::endl;
        summaryFile << "  run time = " << std::fixed << std::setprecision(2)
            << parameters::runTime.count() << " [sec] = "
            << std::fixed << std::setprecision(2) << parameters::runTime.count() / 60 << " [min]" << std::endl;
        summaryFile.close();
    }
    else {
        std::cerr << "Error opening summary file!" << std::endl;
    }
}

void model::Run()
{
    using namespace variables;
    using namespace io;

    Prepare();

    // ----- set the initial conditions of the model -----

    parameters::t = 0.0;
    parameters::step = 0;
    bool solutionFound = true;
    r_au = reader.GetDouble("r_0");
    v_r = reader.GetDouble("v_0");

    // Initial value of the specific angular momentum
    // - the Keplerian angular momentum at the initial point
    j = disk::V_k(r_au) * r_au * astro_const::au / phys_model::scales::Get("j");

    a_nd = reader.GetDouble("a_0") / phys_model::scales::Get("a");
#ifdef STEADY_DRIFT
    sol[0] = r_au;
    sol[1] = a_nd;
#else
    sol[0] = r_au;
    sol[1] = v_r;
    sol[2] = j;
    sol[3] = a_nd;
#endif

    parameters::Print();

    std::cout << "Press Enter to start a simulation" << std::endl;
    std::cin.get();

    // ----- save initial state --------------------------
    SaveCurrentData();

    auto start = std::chrono::high_resolution_clock::now();

    // ----- Main loop ---------------------------------------------
    //       The run stops either when t_stop is reached,
    //       or the grain reaches the inner boundary of the disk
    do
    {
        // ----- read IC by the solver -----
        scheme::solver.ReadInitialVector(sol, parameters::varN);

        // ----- solve the equations -----
        solutionFound = scheme::solver.Solve();
        if (!solutionFound)
		{
			printf("Warning! Solution was not found. Exiting.\n");
			break;
		}

        // ----- update the time step -----

        parameters::dt = scheme::solver.GetTimeStep();
		scheme::solver.UpdateStep(parameters::dt);
        parameters::t += parameters::dt;
		scheme::solver.SetTime(parameters::t);

        // ----- get the solution from the solver -----
		r_au = scheme::solver.GetSolution(0);
        sol[0] = r_au;
#ifdef STEADY_DRIFT
        a_nd = scheme::solver.GetSolution(1);
        sol[1] = a_nd;
        v_r = phys_model::equations::TerminalDriftSpeed(r_au, a_nd * phys_model::scales::Get("a"), disk::GetBeta_k());
        j = disk::V_k(r_au) * r_au * astro_const::au / phys_model::scales::Get("j");
#else
        v_r = scheme::solver.GetSolution(1);
        sol[1] = v_r;
        j = scheme::solver.GetSolution(2);
        sol[2] = j;
        a_nd = scheme::solver.GetSolution(3);
        sol[3] = a_nd;
#endif

		/// save current state
        if ((parameters::step % parameters::savePeriod) == 0)
        {
            std::cout << "Step No " << parameters::step;
            SaveCurrentData();
            std::cout << " - data saved!" << std::endl;
        }

        /// print current step info on screen
        if ((parameters::step % parameters::outputPeriod) == 0)
        {
            std::cout << "Step No " << parameters::step;
            std::cout << std::scientific;
            std::cout << std::setprecision(2);
            std::cout << "\n t  = " << parameters::t << " [t0]" << std::endl;
            std::cout << " dt = " << scheme::solver.GetTimeStep() << " [t0] = " << scheme::solver.GetTimeStep() / steady_solution::tauStopMin << " [tauStopMin]" << std::endl;
            std::cout << " r  = " << r_au << " [au]" << std::endl;
        }

        parameters::step++;

    } while((parameters::t < parameters::stopTime) && (r_au >= io::reader.GetDouble("r_in")));

    auto end = std::chrono::high_resolution_clock::now();
    parameters::runTime = end - start;
    std::cout << "Run time = " << parameters::runTime.count() << " [s]" << std::endl;
    std::cout << "Total number of steps = " << parameters::step << std::endl;

    Finalize();
}
