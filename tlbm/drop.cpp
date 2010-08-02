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
	// Allocate lattice
	LBM::Lattice l(/*FileKey*/"drop", /*Is3D*/false,1.0/6.0, /*Nx*/100, /*Ny*/100,1,1,1);

	// Set constants
	l.SetG(-120.0);
    //l.SetTau(1.0);

	// Set Gravity
	l.SetGravity (0.0, -0.0005);

	// Set walls (top and bottom)
	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid(true, -200);
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid(true, -200);

	// Set drop
	int obsX   = l.Nx()/2; // x position
	int obsY   = l.Ny()/2; // y position
	int radius = l.Nx()/8; // Inital drop radius
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		Vec3_t v0;  v0 = 0.0, 0.0, 0.0;
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
			l.GetCell(i,j)->Initialize (/*Rho*/ 544.0, v0,l.Cs());
		else
			l.GetCell(i,j)->Initialize (/*Rho*/ 8.4e1, v0,l.Cs());
	}

	// Solve
	l.Solve(/*tIni*/0.0, /*tFin*/3000.0, /*dtOut*/10.0);
	return 0;
}
MECHSYS_CATCH
