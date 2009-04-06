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
#include "util/numstreams.h"

using std::cout;
using std::endl;
using Util::_3;
using Util::_6_4;

int main(int argc, char **argv) try
{
	// Allocate lattice
	LBM::Lattice l("basic", /*Is3D*/false,1, /*Nx*/150, /*Ny*/150,1,1,1);
	l.SetG(-0.15);

	// Set output cells
	Array<size_t> outcells;
	outcells.Push(75+75*150);
	l.SetOutCells (outcells,"basic");

	// Set walls (top and bottom)
	//for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	//for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

	// Define Initial conditions: velocity speed and density
	double rho0 = 1.0;
	Vec3_t v0;      v0     = 0.005, 0.0, 0.0;
	Vec3_t vnoise;  vnoise = 0.500, 0.0, 0.0;
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		if (i==75 && j==75) l.GetCell(i,j)->Initialize (rho0, vnoise,l.Cs());
		else                l.GetCell(i,j)->Initialize (rho0, v0,l.Cs());
	}

	// Solve
	l.Solve(/*tIni*/0.0, /*tFin*/200.0, /*dtOut*/1.0);
	//l.Solve(/*tIni*/0.0, /*tFin*/1.0, /*dt*/1.0, /*dtOut*/1.0);
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
