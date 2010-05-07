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
#include <fstream>

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/lbm/lattice.h>

using std::cout;
using std::endl;
using std::ifstream;

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
    if (argc!=2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());
    int         nx;     // Number of cells in the x direction
    int         ny;     // Number of cells in the y direction
    double      Gs;     // Interaction constant for fluid-solid
    double      Gf;     // Interaction constant for fluid-fluid
    double      densv;  // Density for the vapour phase
    double      densl;  // Density for the liquid phase
     
    // Reading from the file
    infile >> nx;     infile.ignore(200,'\n'); 
    infile >> ny;     infile.ignore(200,'\n');
    infile >> Gs;     infile.ignore(200,'\n');
    infile >> Gf;     infile.ignore(200,'\n');
    infile >> densv;  infile.ignore(200,'\n');
    infile >> densl;  infile.ignore(200,'\n');


	// Allocate lattice
	LBM::Lattice l(filekey.CStr(),   // FileKey
	               false,            // Is3D
		           1./6.,
	               nx,               // Nx
	               ny,
		           1,
		           1,
		           1);      

	// Set walls (top and bottom)
	l.SetG(-Gf)->SetGSolid(-Gs);
	l.SetTau(1.0);


	// Define Initial conditions: velocity speed and density
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
	    l.GetCell(i,j)->Initialize (densv, Vec3_t(0.0,0.0,0.0),l.Cs());
	}

	// Set grains
	Table grains;
    String outname (filekey+".out");
	grains.Read(outname.CStr());
	for (size_t i=0; i<grains["Xc"].Size(); ++i)
	{
		double xc = grains["Xc"][i]*nx;
		double yc = grains["Yc"][i]*ny;
		double r  = grains["R" ][i]*0.9*nx;
        if ((xc+r<nx)&&(xc-r>0)&&(yc+r<ny)&&(yc-r>0))
        {
            DrawOpenCircle (l, xc, yc, r  ,0.4);
            DrawFluidCircle(l, xc, yc, r-1,densl);
        }
	}
	l.Solve(/*tIni*/0.0, /*tFin*/40000.0, /*dtOut*/200.);
}
MECHSYS_CATCH
