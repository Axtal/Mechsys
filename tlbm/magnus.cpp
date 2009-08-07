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

// MechSys
#include "lbm/lattice.h"

using std::cout;
using std::endl;

// Analysis constants
double tau    = 1.0;
int    nx     = 100;
int    ny     = 50;
double v0     = 0.04;
int    radius = ny/10;
double omega  = 0.01;  // Angular velocity in the surface of the obstacle

int main(int argc, char **argv) try
{
	// Allocate lattice
	LBM::Lattice l("magnus",      // FileKey
	               false,         // Is3D
				   1,
	               nx,            // Nx
	               ny,
				   1,
				   1,
				   1);           // Ny

	// Set walls (top and bottom)
	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

	// Set inner obstacle
	int obsX = nx/2;   // x position
	int obsY = ny/2; // y position
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		// Calculate the magnus velocity in the obstacle
		double vx = -omega*((int)j-obsY);
		double vy =  omega*((int)i-obsX);

		// Set solids and set the magnus velocity
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
			l.GetCell(i,j)->SetSolid(vx, vy);
	}
	
	// Define boundary conditions
	for (size_t j=0; j<l.Ny(); j++)
	{
		Vec3_t v;  v = v0, 0.0, 0.0;
		l.SetVelocityBC (0, j, v);
		l.SetDensityBC  (nx-1,j, 1.0);
	}

	// Define Initial conditions: velocity speed and density
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		double rho0 = 1.0;
		Vec3_t v;  v = v0, 0.0, 0.0;
		l.GetCell(i,j)->Initialize (rho0, v,l.Cs());
	}

	// Solve
	l.Solve(/*tIni*/0.0, /*tFin*/10000.0, /*dtOut*/10.0);
	//l.Solve(/*tIni*/0.0, /*tFin*/1.0, /*dt*/1.0, /*dtOut*/1.0);
}
MECHSYS_CATCH
