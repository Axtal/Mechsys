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
#include <iostream>
#include <stdlib.h>

// MechSys
#include <mechsys/lbm/Domain.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    size_t nx = 200;
    size_t ny = 200;
    double nu = 0.1;
    double dx = 1.0;
    double dt = 1.0;
    double Tf = 10000.0;
    Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), dx, dt);
    Dom.Lat.G    = -200.0;
    Dom.Lat.Gs   = -200.0;
	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
	{
		
		double rho0 = (200.0 +(1.0*rand())/RAND_MAX)*dx*dx;
		Vec3_t v0;  v0 = 0.0, 0.0, 0.0;
		Dom.Lat.GetCell(iVec3_t(i,j,0))->Initialize (rho0, v0);
	}

	// Solve
    Dom.Solve(Tf,50.0,NULL,NULL,"bubble");
}
MECHSYS_CATCH
