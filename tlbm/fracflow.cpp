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

double distline(double x1,double y1,double x2,double y2,double x,double y) {
	double t,v2;
	v2=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
	t=((x-x1)*(x2-x1)+(y-y1)*(y2-y1))/v2;
	if (t<=0) return sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
	else if (t>=1) return sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2));
	else return sqrt((x-x1-t*(x2-x1))*(x-x1-t*(x2-x1))+(y-y1-t*(y2-y1))*(y-y1-t*(y2-y1)));
	
}

int main(int argc, char **argv) try
{	

	double width=3.;
	// Allocate lattice
	LBM::Lattice l(/*FileKey*/"fracflow", /*Is3D*/false,1, /*Nx*/400, /*Ny*/400,1,1,1);
	
	// Set constants
	l.SetGSolid(0.0);

	// Set walls (top and bottom)
	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		l.GetCell(i,j)->SetSolid();
	}

	// Set grains
	FileParser::Table grains;
	FileParser fp("percolating.txt");
	fp.ReadTable(grains);
	for (size_t i=0; i<grains["x1"].Size(); ++i)
	{
		double x1 = grains["x1"][i]*20.0;
		double y1 = grains["y1"][i]*20.0;
		double x2 = grains["x2"][i]*20.0;
		double y2 = grains["y2"][i]*20.0;
		for (size_t i=0; i<l.Nx(); ++i)
		for (size_t j=0; j<l.Ny(); ++j)
		{
			double x=i;
			double y=j;
			if (distline(x1,y1,x2,y2,x,y)<=width) l.GetCell(i,j)->SetSolid(false);
		}
	}

	
	
	// Initial conditions
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		double      r0 = 1.0;
		Vec3_t v0;  v0 = 0.0, 0.0, 0.0;
		l.GetCell(i,j)->Initialize (r0, v0,l.Cs());
	}

	// Boundary conditions
	for (size_t j=0; j<l.Ny(); j++)
	{
		Vec3_t v; v = 0.01, 0.0, 0.0;
		//l.SetVelocityBC (0,       j, v);
		l.SetDensityBC 	(0,j,1.0);
		l.SetDensityBC  (l.Nx()-1,j, 0.0);
	}
	l.SetTau(1.);
	l.WriteState(0);
	// Solve
	l.Solve(/*tIni*/0.0, /*tFin*/5000.0, /*dtOut*/20.0);
	//l.Solve(/*tIni*/0.0, /*tFin*/1.0, /*dt*/1.0, /*dtOut*/1.0);
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
