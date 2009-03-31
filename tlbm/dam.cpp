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
	LBM::Lattice l(/*FileKey*/"dam", /*Is3D*/false, /*Nx*/150, /*Ny*/150);

	// Set constants
	l.SetTau(0.99)->SetG(-5.5)->SetGSolid(-0.5);

	// Set gravity
	l.SetGravity(0.0,-0.001);

	// Set walls (top and bottom)
	for (size_t i=0; i<l.Top()   .Size(); ++i) l.Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();
	for (size_t i=0; i<l.Right() .Size(); ++i) l.Right()[i]->SetSolid();

	// Set dam and inner obstacle
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		Vec3_t v0;  v0 = 0.0, 0.0, 0.0;
		     if ((i>l.Nx()/2-2)&&(i<(l.Nx()/2+2))&&(j<l.Ny()/10)) l.GetCell(i,j)->SetSolid();
		else if ((i<l.Nx()/3)&&(i>=1)&&(j>=1))                    l.GetCell(i,j)->Initialize(2.6, v0);
		else                                                      l.GetCell(i,j)->Initialize(0.08,v0);
	}
	
	// Solve
	l.WriteState (0);
	l.Solve (/*tIni*/0.0, /*tFin*/1500.0, /*dt*/1.0, /*dtOut*/10.0);
	return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
