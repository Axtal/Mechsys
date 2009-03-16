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
	if (argc>=2) seed = atof(argv[1]);
	srand(seed);

	// Allocate lattice
	LBM::Lattice l(/*FileKey*/ "capi", /*Is3D*/false, /*Tau*/1.0, /*dL*/1.0, /*Nx*/100, /*Ny*/100);

	// Set walls (top and bottom)
	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

	// Set inner walls
	Array<size_t> cols;
	cols.Push(20);
	cols.Push(21);
	cols.Push(30);
	cols.Push(31);
	for (size_t n=0; n<cols.Size(); n++)
	for (size_t j=40; j<60; ++j)
	{
		l.GetCell(cols[n],j)->SetSolid();
	}

	// Initialize cells
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		double rho0 = 1.4 +(.02*rand())/RAND_MAX;
		Vec3_t v0;  v0 = 0.0, 0.0, 0.0;
		l.GetCell(i,j)->Initialize (/*Tau*/1.0, rho0, /*Vx*/v0(0), /*Vy*/v0(1));
	}

	// Solve
	l.SetGravity(0.0, -0.0005);
	l.Solve(/*tIni*/0.0, /*tFin*/400.0, /*dt*/1.0, /*dtOut*/10.0);
	//l.Solve(/*tIni*/0.0, /*tFin*/1.0, /*dt*/1.0, /*dtOut*/1.0);
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
