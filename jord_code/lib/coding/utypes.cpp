#include "utypes.h"

void AllocateVector1D(Tfloat_vector *v, Tint N)
{
	(*v).resize(N);
};

void AllocateVector2D(Tfloat_vector_2d *v, Tint Nx, Tint Ny)
{
	(*v).resize(Nx);
	for (Tint ix = 0; ix < Nx; ix++)
	{
		(*v)[ix].resize(Ny);
	};
};
