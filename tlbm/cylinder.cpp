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
double Re     = 100;
int    nx     = 200;
int    ny     = 50;
int    radius = ny/10 + 1;

double CalcViscosity()
{
	double kin_visc = u_max*(2*radius)/Re; // nu
	return 3.0*kin_visc + 0.5;
}

void CalcInitSpeed(int x, int y, double & vx, double & vy)
{
	double L = ny - 2;
	double yp = y - 1.5;
	vx = u_max*4/(L*L)*(L*yp - yp*yp);
	vy = 0.0;
}

int main(int argc, char **argv) try
{
	// Allocate lattice
	LBM::Lattice l("cylinder",      // FileKey
	               false,           // Is3D
	               nx,              // Nx
	               ny);             // Ny

	// Set tau
	l.SetTau(CalcViscosity());

	// Set walls (top and bottom)
	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

	// Set inner obstacle
	int obsX = ny/2;   // x position
	int obsY = ny/2+3; // y position
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
			l.GetCell(i,j)->SetSolid();
	}
	
	// Define boundary conditions
	for (size_t j=0; j<l.Ny(); j++)
	{
		double vx, vy;
		CalcInitSpeed (0, j, vx, vy);
		Vec3_t v;  v = vx, vy, 0.0;
		l.SetVelocityBC (0, j, v);
		l.SetDensityBC  (nx-1,j, 1.0);
	}

	// Define Initial conditions: velocity speed and density
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		double rho0 = 1.0;
		Vec3_t v0;  v0 = 0.08, 0.0, 0.0;
		l.GetCell(i,j)->Initialize (rho0, v0);
	}

	// Solve
	l.Solve(/*tIni*/0.0, /*tFin*/15000.0, /*dt*/1.0, /*dtOut*/200.0);
	//l.Solve(/*tIni*/0.0, /*tFin*/1.0, /*dt*/1.0, /*dtOut*/1.0);
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
