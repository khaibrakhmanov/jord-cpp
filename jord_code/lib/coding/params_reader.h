#pragma once

#include <string>
#include <map>

/// @brief Class - reader of parameters from a text file
class TParameterReader {
public:
    /// @brief Read the parameters file
    /// @param fname - file name
    void ReadParameters(const std::string& fname);

    /// @brief Get the double-type parameter from the parameters dictionary
    /// @param name - parameter name
    /// @return parameter value
    double GetDouble(const std::string& name) const;

    /// @brief Get the int-type parameter from the parameters dictionary
    /// @param name - parameter name
    /// @return parameter value
    int GetInt(const std::string& name) const;

private:
    /// @brief the dictionary containing the parameters read
    std::map<std::string, std::string> paramsDict;
};
