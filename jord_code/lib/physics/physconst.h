#pragma once

#include "../coding/utypes.h"


/* fundamental physical constants */
namespace physical_const
{
	/* в единицах СГС */

	/// 1 год, с
	const Tfloat year = 3.1536e7;
	/// скорость света, см/с
	const Tfloat c = 2.99792458e10;
	/// постоянная Больцмана
	const Tfloat k = 1.38064852e-16;
	/// заряд электрона
	const Tfloat e = 4.803204673e-10;
	/// универсальная газовая постоянная
	const Tfloat R = 8.3144598e7;
	/// гравитационная постоянная
	const Tfloat G = 6.674e-8;
	/// боровский радиус атома водорода
	const Tfloat aH = 5.29e-9;
	/// масса протона
	const Tfloat m_p = 1.672621e-24;
	/// постоянная Стефана-Больцмана
	const Tfloat sigma_R = 5.670367e-5;
	/// масса электрона
	const Tfloat m_e = 9.1e-28;
	/// электронвольт, эрг
	const Tfloat  eV = 4.8e-10;
	/// постоянная Планка, эрг*с
	const Tfloat  h = 6.62e-27;;

	/* in SI units */
	namespace SI
	{
		/* speed of light */
		const Tfloat c = 2.99792458e8;
	}

}