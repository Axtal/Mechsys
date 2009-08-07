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
int    nx     = 200;
int    ny     = 200;
double v0     = 0;
//double omega  = 0.001;  // Angular velocity in the surface of the obstacle


void DrawCircle(LBM::Lattice & l, double & obsX, double & obsY, double radius, double vx, double vy, double & fx, double & fy, double dt, int T)
{
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		// Set solids and set the magnus velocity
		if (pow((int)obsX-(int)i,2.0) + pow((int)obsY-(int)j,2.0) <= pow(radius,2.0)) // circle equation
		{
			l.GetCell(i,j)->SetSolid(vx, vy);
		}
		//else 
			//l.GetCell(i,j)->SetSolid(false);
	}

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
				//size_t op = c->Opp(k);
				//double alpha = 6.0*c->W(op)*rho;
				//fx += (2.0*c->F(op) - alpha*(c->C(op,0)*vx+c->C(op,1)*vy ))/dt*c->C(op,0);
				//fy += (2.0*c->F(op) - alpha*(c->C(op,0)*vx+c->C(op,1)*vy ))/dt*c->C(op,1);

				double G = -3.0;
				double rho_nb = l.GetCell(c->Neigh(k))->Density();
				double nb_psi = l.Psi(rho_nb);
				fx += -G*c->W(k)*nb_psi*c->C(k,0);
				fy += -G*c->W(k)*nb_psi*c->C(k,1);
			}
		}
	}

}

int main(int argc, char **argv) try
{
	// Allocate lattice
	LBM::Lattice l("cap_grain",   // FileKey
	               false,    // Is3D
		       1./6.,
	               nx,       // Nx
	               ny,
		       1,
		       1,
		       1);      

	// Set walls (top and bottom)
	l.SetG(-6.0)->SetGSolid(-3.5);
	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

	//// Define boundary conditions
	/*for (size_t j=0; j<l.Ny(); j++)
	{
		Vec3_t v;  v = 0.0, 0.0, 0.0;
		l.SetVelocityBC (0, j, v);
		l.SetDensityBC  (nx-1,j, 1.0);
	}
	*/
	// Define Initial conditions: velocity speed and density
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		//double rho0  = 1.0;
		size_t radio = 24;
		Vec3_t V;  V = 0.0, 0.0, 0.0;
		//if (pow((int)(i)-nx/2,2.0) + pow((int)(j)-20,2.0) <= pow(radio,2.0)) // circle equation
		if (pow((int)(i)-nx/2,2.0) + pow((int)(j)-100,2.0) <= pow(radio,2.0)) // circle equation
		{
			l.GetCell(i,j)->Initialize (1.5, V,l.Cs());
		}
		else
		{
			l.GetCell(i,j)->Initialize (0.05, V,l.Cs());
		}
	}

	size_t Tmax = 2000;
	size_t Tout = 10;

	double vx = 0.0;
	double vy = 0.0;
	double dt = 1.0;
	double m  = 400.;
	// Set inner obstacle
	double obsX   = nx/2; // x position
	double obsY   = ny/3; // y position
	int    radius = 30; //ny/10;
	double fx     = 0.0;
	double fy     = 0.0;
	double X1,Y1,X2,Y2,vx1,vy1,vx2,vy2,fx1,fy1,fx2,fy2;
	X1=nx/3;
	X2=2*nx/3;
	Y1=Y2=2*ny/3;
	vx1=vx2=vy1=vy2=fx1=fx2=fy1=fy2=0.0;

	DrawCircle(l, obsX, obsY, radius, vx, vy, fx, fy, dt,0);
	DrawCircle(l, X1, Y1,radius,vx1,vy1,fx1,fy1,dt,0);
	DrawCircle(l, X2, Y2,radius,vx2,vy2,fx2,fy2,dt,0);

	for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
	for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

	std::cout << "Stage 0" << std::endl;
	l.WriteState (0);

	for (size_t T=1; T<=Tmax; T++)
	{
		for (size_t i=0; i<l.Nx(); i++)
			for (size_t j=0; j<l.Ny(); j++)
				{
					l.GetCell(i,j)->SetSolid(false);
				}

		fx = fx1=fx2=0.0;
		fy = fy1=fy2=0.0;

		DrawCircle(l, obsX, obsY, radius, vx, vy, fx, fy, dt,T);
		DrawCircle(l, X1, Y1,radius,vx1,vy1,fx1,fy1,dt,T);
		DrawCircle(l, X2, Y2,radius,vx2,vy2,fx2,fy2,dt,T);
		for (size_t i=0; i<l.Top()   .Size(); ++i) l   .Top()[i]->SetSolid();
		for (size_t i=0; i<l.Bottom().Size(); ++i) l.Bottom()[i]->SetSolid();

		if (T>1000)
		{
			vx += fx/m*dt;
			vy += fy/m*dt;
			vx1 += fx1/m*dt;
			vy1 += fy1/m*dt;
			vx2 += fx2/m*dt;
			vy2 += fy2/m*dt;
			obsX += vx*dt;
			obsY += vy*dt;
			X1 += vx1*dt;
			X2 += vx2*dt;
			Y1 += vy1*dt;
			Y2 += vy2*dt;
		}

		//std::cout << "obsX = " << obsX << std::endl;
		//std::cout << "obsY = " << obsY << std::endl;
		//std::cout << "fx = " << fx << std::endl;
		//std::cout << "fy = " << fy << std::endl;
		//std::cout << "vx = " << vx << std::endl;
		//std::cout << "vy = " << vy << std::endl;

		l.ApplyForce   ();
		l.Collide      ();
		l.BounceBack   ();
		l.Stream       ();
		l.ApplyBC      ();
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
