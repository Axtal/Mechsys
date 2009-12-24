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
#include <mechsys/lbm/lattice.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	// Input 
	double seed = 100.0;
	double h=1.,dt=1.;
	if (argc>=2) seed = atof(argv[1]);
	srand(seed);

	// Allocate lattice
	LBM::Lattice l("bubble", // FileKey
	               false,    // Is3D
	               1./6., 	 //viscosity
	               int(75/h),// Nx
	               int(75/h),// Ny
		       1, 	 // Nz
		       h,	 // h
		       dt   	 // dt
		       );     

	// Set constants
	std::cout << l.Tau() << " "<< l.dt()<<std::endl;
	l.SetG(-6.0);

	// Initialize cells
	
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		
		double rho0 = (0.5 +(.1*rand())/RAND_MAX)*h*h;
		Vec3_t v0;  v0 = 0.0, 0.0, 0.0;
		l.GetCell(i,j)->Initialize (rho0, v0,l.Cs());
		//std::cout<<l.GetCell(i,j)->Density()<<std::endl;

	}
	
	
	/*for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		double rho0;
		if ((i-l.Nx()/2)*(i-l.Nx()/2)+(j-l.Ny()/2)*(j-l.Ny()/2)<=400) 
			rho0 = 2.;
		else rho0=1.;
		Vec3_t v0;  v0 = 0.0, 0.0, 0.0;
		l.GetCell(i,j)->Initialize (rho0, v0,l.Cs());
		//std::cout<<l.GetCell(i,j)->Density()<<std::endl;

	}
	*/

	// Solve
	l.Solve(/*tIni*/0.0, /*tFin*/2000.0, /*dtOut*/1.);
}
MECHSYS_CATCH
