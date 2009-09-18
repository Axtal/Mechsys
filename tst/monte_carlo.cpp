/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// Std Lib
#include <stdlib.h>
#include <math.h>

// MechSys
#include "util/fatal.h"
#include "numerical/montecarlo.h"

using Numerical::MonteCarlo;

class sphere
{
	public:
		double Inside(double * r)
		{
			if ((r[0]-x)*(r[0]-x)+(r[1]-y)*(r[1]-y)+(r[2]-z)*(r[2]-z)<R*R) return 1.;
			else return 0.;
		}
		double dInertia(double *r)
		{
			return (r[0]*r[0]+r[1]*r[1])*Inside(r);
		}
		double x,y,z,R;
};

int main(int argc, char **argv) try
{
	double ri[3] = { -5, -5, -5 };
	double rs[3] = { 5, 5, 5 };
	double sr = 3;
	sphere S;
	S.x = S.y = S.z = 0;
	S.R = sr;
	MonteCarlo<sphere> MC1;
	std::cout << " Exact Volume = "<<(4./3.)*M_PI*sr*sr*sr << " Calculated Volume = " <<MC1.Integrate(&S,&sphere::Inside,ri,rs,Numerical::VEGAS,500000)<<std::endl;
	std::cout << " Exact Inertia Moment = "<<(8./15.)*M_PI*sr*sr*sr*sr*sr << " Calculated Inertia moment = " <<MC1.Integrate(&S,&sphere::dInertia,ri,rs,Numerical::VEGAS,500000)<<std::endl;

	return 0;
}
MECHSYS_CATCH
