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
#include "util/fileparser.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	// Allocate lattice
	LBM::Lattice l(/*FileKey*/"porous", /*Is3D*/false, /*Nx*/200, /*Ny*/200);
	
	// Set constants
	l.SetTau(1.0)->SetGSolid(0.0);

	// Set walls (top and bottom)
	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

	// Set grains
	FileParser::Table grains;
	FileParser fp("circles.out");
	fp.ReadTable(grains);
	for (size_t i=0; i<grains["Xc"].Size(); ++i)
	{
		double xc = grains["Xc"][i]*200.0;
		double yc = grains["Yc"][i]*200.0;
		double r  = grains["R" ][i]*190.0;
		for (size_t i=0; i<l.Nx(); ++i)
		for (size_t j=0; j<l.Ny(); ++j)
		{
			double dd = pow(static_cast<double>(i)-xc,2.0) + pow(static_cast<double>(j)-yc,2.0);
			if (dd<=r*r) l.GetCell(i,j)->SetSolid();
		}
	}
	
	// Boundary conditions
	for (size_t j=0; j<l.Ny(); j++)
	{
		Vec3_t v; v = 0.01, 0.0, 0.0;
		l.SetVelocityBC (0,       j, v);
		l.SetDensityBC  (l.Nx()-1,j, 1.0);
	}

	// Initial conditions
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		double      r0 = 1.0;
		Vec3_t v0;  v0 = 0.001, 0.0, 0.0;
		l.GetCell(i,j)->Initialize (r0, v0);
	}

	// Solve
	l.Solve(/*tIni*/0.0, /*tFin*/15000.0, /*dt*/1.0, /*dtOut*/20.0);
	//l.Solve(/*tIni*/0.0, /*tFin*/1.0, /*dt*/1.0, /*dtOut*/1.0);
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
