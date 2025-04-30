#ifndef ASTROCONST_INCLUDED
#define ASTROCONST_INCLUDED

#include "../coding/utypes.h"

/*!
  \file astroconst.h
  \brief Астрономические константы
*/

/*!
  \namespace astro_const
  \brief Пространство имен для астрономических констант 
*/
namespace astro_const
{
	/// астрономическая единица, см
	const Tfloat au = 1.5e13;
	/// радиус Солнца, см
	const Tfloat R_sun = 6.96e10;
	/// масса Солнца, г
	const Tfloat M_sun = 1.99e33;
	/// средний молекулярный вес газа в МЗС
	const Tfloat mu_ISM = 2.3;
	///год, с
	const Tfloat year = 365.0 * 24.0 * 3600.0;
};

#endif