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

void DrawBullEyeCircle(LBM::Lattice & l, double obsX, double obsY, double radius, double thickness, size_t n_div, double arc)
{
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
        double r2 = pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0);
		if ((r2>(1-thickness)*(1-thickness)*radius*radius)&&(r2<radius*radius))
		{
            double angle = atan2(int(j)-obsY,int(i)-obsX)*180.0/M_PI;
            bool solid = true;
            for (size_t k=0;k<n_div;k++)
            {
                double angled = -180.0 + k*360.0/n_div;
                if (fabs(angle-angled)<arc) solid = false;
            }
			if (solid) l.GetCell(i,j)->SetSolid();
		}
    }
}

void DrawFluidDisk(LBM::Lattice & l , double obsX, double obsY, double radius, double density)
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

void DrawFluidCircle(LBM::Lattice & l , double obsX, double obsY, double r1, double r2, double density)
{
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
        double r = pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0);
		if (r>=r1*r1&&r<r2*r2) // circle equation
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
    int         Tf;     // final time
    int         dTout;  // time span for output
    int         n_div;  // number of divitions
    double      arc;    // arc length of holes
     
    // Reading from the file
    infile >> nx;     infile.ignore(200,'\n'); 
    infile >> ny;     infile.ignore(200,'\n');
    infile >> Gs;     infile.ignore(200,'\n');
    infile >> Gf;     infile.ignore(200,'\n');
    infile >> densv;  infile.ignore(200,'\n');
    infile >> densl;  infile.ignore(200,'\n');
    infile >> Tf;     infile.ignore(200,'\n');
    infile >> dTout;  infile.ignore(200,'\n');
    infile >> n_div;  infile.ignore(200,'\n');
    infile >> arc;    infile.ignore(200,'\n');


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
            //DrawOpenCircle (l, xc, yc, r+1,0.4);
            DrawBullEyeCircle (l, xc, yc, 1.1*r,0.1,n_div,arc);
            DrawFluidDisk(l, xc, yc, 0.9*1.1*r ,densl);
            DrawFluidCircle(l, xc, yc, 1.1*r,1.3*r, densl);
        }
	}
	l.Solve(/*tIni*/0.0, /*tFin*/Tf, /*dtOut*/dTout);
}
MECHSYS_CATCH
