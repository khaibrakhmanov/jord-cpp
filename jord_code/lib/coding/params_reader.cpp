#include "params_reader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

void TParameterReader::ReadParameters(const std::string& fname) {
    std::ifstream inputFile(fname);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Could not open parameter file: " + fname);
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        // Skip comments and empty lines
        if (line.substr(0, 2) == "//" || line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        std::string type, name, eq, value;
        if (iss >> type >> name >> eq >> value && eq == "=") {
            paramsDict[name] = value;
        } else {
            throw std::runtime_error("Invalid line format in parameter file: " + line);
        }
    }

    inputFile.close();
}

double TParameterReader::GetDouble(const std::string& name) const {
    auto it = paramsDict.find(name);
    if (it == paramsDict.end()) {
        throw std::runtime_error("Parameter not found: " + name);
    }
    std::istringstream iss(it->second);
    double value;
    iss >> value;
    if (iss.fail()) {
        throw std::runtime_error("Could not convert parameter '" + name + "' to the requested type.");
    }
    return double(value);
}

int TParameterReader::GetInt(const std::string& name) const {
    auto it = paramsDict.find(name);
    if (it == paramsDict.end()) {
        throw std::runtime_error("Parameter not found: " + name);
    }
    std::istringstream iss(it->second);
    int value;
    iss >> value;
    if (iss.fail()) {
        throw std::runtime_error("Could not convert parameter '" + name + "' to the requested type.");
    }
    return int(value);
}
