#include <iostream>
#include "num_model.h"

using namespace std;

int main()
{
    cout << "Hello! The program performs the numerical simulation of the dust grain dynamics." << endl;

    // read model parameters
    model::ReadParameters("params.ini");
    // run the simulation
    model::Run();

    cout << "Simulation completed" << endl;

    return 0;
}
