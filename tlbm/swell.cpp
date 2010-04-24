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
#include <mechsys/lbm/dem.h>

using std::cout;
using std::endl;
int    nx     = 200;
int    ny     = 200;
//double omega  = 0.001;  // Angular velocity in the surface of the obstacle


void DrawOpenCircle(LBM::Lattice & l, double obsX, double obsY, double radius, double tol)
{
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		if (fabs(sqrt(pow((int)obsX-(int)i,2.0) + pow((int)obsY-(int)j,2.0))-radius)<tol) // circle equation
		{
			l.GetCell(i,j)->SetSolid();
		}
	}
}

void DrawFluidCircle(LBM::Lattice & l , double obsX, double obsY, double radius, double density)
{
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
		{
			l.GetCell(i,j)->Initialize (density, Vec3_t(0.0,0.0,0.0),l.Cs());
		}
    }
}

int main(int argc, char **argv) try
{
	// Allocate lattice
	LBM::Lattice l("swell",   // FileKey
	               false,    // Is3D
		           1./6.,
	               nx,       // Nx
	               ny,
		           1,
		           1,
		           1);      

	// Set walls (top and bottom)
	l.SetG(-6.0)->SetGSolid(-10.5);

	// Define Initial conditions: velocity speed and density
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
	    l.GetCell(i,j)->Initialize (0.1, Vec3_t(0.0,0.0,0.0),l.Cs());
	}
    DrawOpenCircle (l,100, 80,11,0.2);
    DrawFluidCircle(l,100, 80,10,2.8);

    DrawOpenCircle (l,120, 90,11,0.2);
    DrawFluidCircle(l,120, 90,10,2.8);

    DrawOpenCircle (l, 80,110,11,0.2);
    DrawFluidCircle(l, 80,110,10,2.8);

    DrawOpenCircle (l,100,120,11,0.2);
    DrawFluidCircle(l,100,120,10,2.8);


	l.Solve(/*tIni*/0.0, /*tFin*/5000.0, /*dtOut*/50.);
}
MECHSYS_CATCH
