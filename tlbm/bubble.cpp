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
#include "lbm/lattice.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	// Input 
	double seed = 5.0;
	double h=1.,dt=0.5;
	if (argc>=2) seed = atof(argv[1]);
	srand(seed);

	// Allocate lattice
	LBM::Lattice l("bubble", // FileKey
	               false,    // Is3D
	               0.1, 	 //viscosity
	               int(75/h),// Nx
	               int(75/h),// Ny
		       1, 	 // Nz
		       h,	 // h
		       dt   	 // dt
		       );     

	// Set constants
	std::cout << l.Tau() << " "<< l.dt()<<std::endl;
	l.SetG(-5.0);

	// Initialize cells
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		double rho0 = (0.8 +(.1*rand())/RAND_MAX)*h*h;
		Vec3_t v0;  v0 = 0.0, 0.0, 0.0;
		l.GetCell(i,j)->Initialize (rho0, v0,l.Cs());
	}

	// Solve
	l.Solve(/*tIni*/0.0, /*tFin*/5000.0, /*dtOut*/1.);
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
