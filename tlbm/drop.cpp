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
	// Allocate lattice
	LBM::Lattice l(/*FileKey*/ "drop", /*Is3D*/false, /*Tau*/1.2, /*dL*/1.0, /*Nx*/100, /*Ny*/100);

	srand(5.0);

	// Set walls (top and bottom)
	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

	// Set inner obstacle
	int obsX = l.Nx()/2;   // x position
	int obsY = l.Ny()/4;   // y position
	int radius = l.Nx()/8; // Inital drop radius

	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
			l.GetCell(i,j)->Initialize (/*Tau*/1.0,/*Rho*/ 0.8, /*Vx*/0.0, /*Vy*/0.0);
		else
			l.GetCell(i,j)->Initialize (/*Tau*/1.0,/*Rho*/ 0.1, /*Vx*/0.0, /*Vy*/0.0);
	}

	// Solve
	l.Solve(/*tIni*/0.0, /*tFin*/200.0, /*dt*/1.0, /*dtOut*/10.0);

	l.SetGravity(0.0, -0.0005);

	l.Solve(/*tIni*/0.0, /*tFin*/1500.0, /*dt*/1.0, /*dtOut*/10.0);


}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
