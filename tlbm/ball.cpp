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


void DrawBall(LBM::Lattice & l, double & obsX, double & obsY, double radius, double vx, double vy, double & fx, double & fy, double dt, int T)
{
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		// Set solids and set the magnus velocity
		if (pow((int)obsX-(int)i,2.0) + pow((int)obsY-(int)j,2.0) <= pow(radius,2.0)) // circle equation
			l.GetCell(i,j)->SetSolid(vx, vy);
		else 
			l.GetCell(i,j)->SetSolid(false);
	}

	// Calculate de force in the solid
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		LBM::Cell * c = l.GetCell(i,j);
		if (c->IsSolid())
		{
			double rho = 0.0;
			for (size_t k=0; k<l.NNeigh(); k++) rho += c->F(k);
			
			for (size_t k=0; k<l.NNeigh(); ++k)
			{
				if (l.GetCell(c->Neigh(k))->IsSolid()) continue;
				size_t op = c->Opp(k);
				double alpha = 6.0*c->W(op)*rho;
				fx += (2.0*c->F(op) - alpha*(c->C(op,0)*vx+c->C(op,1)*vy ))/dt*c->C(op,0);
				fy += (2.0*c->F(op) - alpha*(c->C(op,0)*vx+c->C(op,1)*vy ))/dt*c->C(op,1);
			}
		}
	}
}

int main(int argc, char **argv) try
{
	// Analysis constants
	//double tau    = 1.0;
	int    nx     = 50;
	int    ny     = 20;
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

	//l.SetTau(tau);

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

	// Set inner ball
	double obsX   = ny/2;    // center x position
	double obsY   = ny/2+3;  // center y position
	int    radius = 4; 

	// Initial velocity and force
	double vx = 0.0;
	double vy = 0.0;
	double fx = 0.0;
	double fy = 0.0;

	// Time and mass
	double dt = 1.0;     // Time increment (equal to 1.00 for LB models)
	double m  = 100.0;   // Total mass of the ball

	DrawBall(l, obsX, obsY, radius, vx, vy, fx, fy, dt,0);

	std::cout << "Stage 0" << std::endl;
	l.WriteState (0);

	for (size_t T=1; T<=Tmax; T++)
	{
		// Reset the force
		fx = 0.0;
		fy = 0.0;

		DrawBall(l, obsX, obsY, radius, vx, vy, fx, fy, dt,T);

		// Update the velocity and the position of the ball
		vx   += fx/m*dt;
		vy   += fy/m*dt;
		obsX += vx*dt;
		obsY += vy*dt;

		// Screen output
		/*
		std::cout << "obsX = " << obsX << std::endl;
		std::cout << "obsY = " << obsY << std::endl;
		std::cout << "fx   = " << fx << std::endl;
		std::cout << "fy   = " << fy << std::endl;
		std::cout << "vx   = " << vx << std::endl;
		std::cout << "vy   = " << vy << std::endl;
		*/

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
