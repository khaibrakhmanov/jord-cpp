#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

#include <string>

/// @brief Namespace - interface to the numerical model of dust grain dynamics
namespace model
{
    /// @brief Run the simulation
    void Run();

    /// @brief Read model parameters from the ini-file
    /// @param iniFile - file name (with relative path)
    void ReadParameters(std::string iniFile);
}
#endif // MODEL_H_INCLUDED
