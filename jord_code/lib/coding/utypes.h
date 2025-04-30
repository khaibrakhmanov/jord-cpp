#ifndef UTYPES_INCLUDED
#define UTYPES_INCLUDED


#include <vector>
#include <functional>
#include <map>
#include <string>

/*
	unit contains user-defined types
*/


// user type: float number
typedef double Tfloat;

// user type: integer
typedef unsigned int Tint;

// user type: short integer
typedef short Tshort;

// user type: boolean
typedef bool Tbool;

// user type: one-dimensional float vector
typedef std::vector<double> Tfloat_vector;

// user type: two-dimensional float vector
typedef std::vector<Tfloat_vector> Tfloat_vector_2d;

/// @brief Allocate memory for 1D-vector
/// @param v - pointer to the vector
/// @param N - number of elements
void AllocateVector1D(Tfloat_vector *v, Tint N);

/// @brief Allocate memory for 2D-vector
/// @param v - pointer to the vector
/// @param Nx - number of elements in the first dimension
/// @param Ny - number of elements in the second dimension
void AllocateVector2D(Tfloat_vector_2d *v, Tint Nx, Tint Ny);

/// Dictionary: [name <-> floating-point-type value]
typedef std::map <std::string, Tfloat> TFloatValueDict;

#endif