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
using LBM::Disk;

// Analysis constants
double u_max  = 0.1;
double Re     = 100;
int    nx     = 400;
int    ny     = 400;
int    radius = 5;

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

	double h=1;
	double dt=1;

	// Allocate lattice
	LBM::Lattice l("colloid", false, u_max*(2*radius)/Re, nx, ny, 1, h, dt);
	l.SetTau(1.0);



    // set balls
    Array<Disk> Ball;
    for(double y = radius+10; y <=ny;y+=radius+10)
    {
        Ball.Push(Disk(Vec3_t(nx/20, y, 0.0), Vec3_t(0.0,0.0,0.0), radius, 400.0, 10.0, dt));
    }
					

	// Set walls
    Array<Disk> Wall;
    Wall.Push(Disk(Vec3_t((0.45*double(nx)), 480  , 0.0), Vec3_t(0.0,0.0,0.0), 120, 50.0, 10.0, dt));
    Wall.Push(Disk(Vec3_t((0.45*double(nx)), -80  , 0.0), Vec3_t(0.0,0.0,0.0), 120, 50.0, 10.0, dt));
    Wall.Push(Disk(Vec3_t((0.45*double(nx)), 200  , 0.0), Vec3_t(0.0,0.0,0.0), 150, 50.0, 10.0, dt));

    for (size_t j=0;j<Wall.Size();j++)
    {
        Wall[j].DrawDisk(l,dt);
    }

	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();
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
		l.GetCell(i,j)->Initialize (rho0, v0, l.Cs());
	}


	size_t Tmax   = 30000; // max number of steps
	size_t Tout   = 100;    // interval steps for output

	// Solve
	for (size_t T=0; T<=Tmax; T++)
	{
        l.SetSolid(false);

		// Reset the force
        for (size_t i=0;i<Ball.Size();i++)
        {
            Ball[i].StartForce(Vec3_t(-0.1*Ball[i].v));
            Ball[i].DrawDisk(l,dt);
            for (size_t j=0;j<Wall.Size();j++)
            {
                CalcForce(&Ball[i],&Wall[j]);
            }
        }
        for (size_t j=0;j<Wall.Size();j++)
        {
		    Wall[j].DrawDisk(l,dt);
        }

	    for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	    for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

        for (size_t i=0;i<Ball.Size()-1;i++)
        {
            for (size_t j=i+1;j<Ball.Size();j++)
            {
                CalcForce(&Ball[i],&Ball[j]);
            }
        }
        for (size_t i=0;i<Ball.Size();i++)
        {
            Ball[i].Move(dt);
        }

		l.ApplyForce   ();
		l.Collide      ();
		l.BounceBack   ();
		l.Stream       ();
		l.ApplyBC      ();

		// File output
		if (T % Tout == 0) 
		{
			String buf;
			buf.Printf("[1;34mMechSys[0m::LBM::Lattice::Solve: [1;31mt = %d   TotalMass = %g[0m\n", T, l.TotalMass());
			std::cout << buf;
			l.WriteState (T);
		}
	}
}
MECHSYS_CATCH
