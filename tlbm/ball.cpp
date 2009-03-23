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
double tau    = 1.0;
int    nx     = 3;
int    ny     = 1;
double v0     = 0.01;
//double omega  = 0.001;  // Angular velocity in the surface of the obstacle

void DrawBall(int ic, int jc)
{

}

int main(int argc, char **argv) try
{
	// Allocate lattice
	LBM::Lattice l("ball",      // FileKey
	               false,         // Is3D
	               nx,            // Nx
	               ny);           // Ny

	// Set walls (top and bottom)
	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

	// Define boundary conditions
	for (size_t j=0; j<l.Ny(); j++)
	{
		Vec3_t v;  v = v0, 0.0, 0.0;
		l.SetVelocityBC (0, j, v);
		l.SetDensityBC  (nx-1,j, 1.0);
	}

	// Define Initial conditions: velocity speed and density
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		double rho0 = 1.0;
		Vec3_t v;  v = v0, 0.0, 0.0;
		l.GetCell(i,j)->Initialize (rho0, v);
		for (size_t k=0; k<l.NNeigh(); ++k)
			if (k==1&&i==0) l.GetCell(i,j)->F(k)=1;
			else  l.GetCell(i,j)->F(k)=0.001;


				
	}

	size_t Tmax = 2;
	//size_t Tout = 10;

	double vx = 0.0;
	double vy = 0.0;
	double dt = 0.01;
	double m  = 5;
	// Set inner obstacle
	double obsX   = nx/2; // x position
	double obsY   = ny/2; // y position
	int radius = 0; //ny/10;

	l.WriteState (0);
	for (size_t T=0; T<Tmax; T++)
	{
		double fx = 0.0;
		double fy = 0.0;
		for (size_t i=0; i<l.Nx(); ++i)
		for (size_t j=0; j<l.Ny(); ++j)
		{
			// Set solids and set the magnus velocity
			if (pow((int)obsX-(int)i,2.0) + pow((int)obsY-(int)j,2.0) <= pow(radius,2.0)) // circle equation
			{
				LBM::Cell *   c = l.GetCell(i,j);
				double rho = 0.0;
				for (size_t k=0; k<l.NNeigh(); k++) rho += c->F(k);
				
				for (size_t k=0; k<l.NNeigh(); ++k)
				{
					if (l.GetCell(c->Neigh(k))->IsSolid()) continue;
					size_t op = c->Opp(k);
					double alpha = 6.0*c->W(op)*rho;
					fx += 2.0*(c->F(op) - alpha*(c->C(op,0)*vx+c->C(op,1)*vy ))/dt*c->C(op,0);
					fy += 2.0*(c->F(op) - alpha*(c->C(op,0)*vx+c->C(op,1)*vy ))/dt*c->C(op,1);
				}
				l.GetCell(i,j)->SetSolid(vx, vy);
			}
			else 
				l.GetCell(i,j)->SetSolid(false);
		}

		if (T>100)
		{
			vx += fx/m*dt;
			vy += fy/m*dt;

			obsX += vx*dt;
			obsY += vy*dt;
		}

		std::cout << "obsX = " << obsX << std::endl;
		std::cout << "obsY = " << obsY << std::endl;
		std::cout << "fx = " << fx << std::endl;
		std::cout << "fy = " << fy << std::endl;


		for (size_t k=0; k<l.NNeigh(); k++) 
		{
			//if( l.GetCell(obsX-radius, ny/2)->IsSolid())
				std::cout << "F0(" << k << ") = " << l.GetCell(0,0)->F(k) << std::endl;
				std::cout << "F1(" << k << ") = " << l.GetCell(1,0)->F(k) << std::endl;
				std::cout << "F2(" << k << ") = " << l.GetCell(2,0)->F(k) << std::endl;
		}

		l.ApplyForce   ();
		l.ApplyGravity ();
		l.Collide      ();
		l.BounceBack   ();
		l.Stream       ();
		l.ApplyBC      ();
		//if (T % Tout == 0) 
		{
			String buf;
			buf.Printf("[1;34mMechSys[0m::LBM::Lattice::Solve: [1;31mt = %d   TotalMass = %g[0m\n", T, l.TotalMass());
			std::cout << buf;
			l.WriteState (T);
		}
		
	}

}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
