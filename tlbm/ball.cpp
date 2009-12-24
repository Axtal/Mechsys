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


int main(int argc, char **argv) try
{
	// Analysis constants
	double tau    = 1.0;
	int    nx     = 50;
	int    ny     = 50;
	double v0     = 0.1;
	size_t Tmax   = 3000; // max number of steps
	size_t Tout   = 1;    // interval steps for output

	// Allocate lattice
	LBM::Lattice l("ball",      // FileKey
	               false,       // Is3D
				   1.0,			// viscosity
	               nx,          // Nx
	               ny, 			// Ny
				   1,           // Nz
				   1, 			// space step
				   1 			// time step
				   );         

	l.SetTau(tau);

	// Define velocity and density boundary conditions
	for (size_t j=0; j<l.Ny(); j++)
	{
		Vec3_t v;  v = j*v0/ny, 0.0, 0.0;
		l.SetVelocityBC (0, j, v);
		l.SetDensityBC  (nx-1,j, 1.0);
	}

	// Define Initial conditions: velocity speed and density
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		double rho0 = 1.0;
		Vec3_t v;  v = v0, 0.0, 0.0;
		l.GetCell(i,j)->Initialize (rho0, v,l.Cs());
	}



	// Time and mass
	double dt = 1.0;     // Time increment (equal to 1.00 for LB models)

    LBM::Disk Ball(Vec3_t(nx/3.0, ny/2.0, 0.0), Vec3_t(0.0,0.0,0.0), 4.0, 100.0, 1000.0, dt);
    Ball.DrawDisk(l,dt);

	l.WriteState (0);

	for (size_t T=1; T<=Tmax; T++)
	{
		// Reset the force
        Ball.StartForce();
        //l.SetSolid(false);
        Ball.DrawDisk(l,dt);

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
