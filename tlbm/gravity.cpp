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
double u_max  = 0.1;
double Re     = 30;
int    nx     = 50;
int    ny     = 50;
int    radius = ny/10 + 1;

int main(int argc, char **argv) try
{
	// Allocate lattice
	LBM::Lattice l(/*FileKey*/ "grav", /*Is3D*/false, /*Tau*/1.0, /*dL*/1.0, nx, ny);

	// Set walls (top and bottom)
	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

	//// Set inner obstacle
	//int obsX = nx/2;   // x position
	//int obsY = ny/2+3; // y position
	//for (size_t i=0; i<l.Nx(); ++i)
	//for (size_t j=0; j<l.Ny(); ++j)
	//{
	//	if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
	//		l.GetCell(i,j)->SetSolid();
	//}
	
	// Define Initial conditions: velocity speed and density
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		double tau = 1.0;
		double rho = 1.0;
		//l.GetCell(i,j)->Initialize (tau, rho, vx, vy);
		l.GetCell(i,j)->Initialize (tau, rho, 0.0, 0.0);
	}

	// Solve
	l.Solve(/*tIni*/0.0, /*tFin*/10000.0, /*dt*/1.0, /*dtOut*/20.0);
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }