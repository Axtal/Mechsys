/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#ifndef MECHSYS_LBM_DEM_H
#define MECHSYS_LBM_DEM_H

#include "lbm/mixture.h"

namespace LBM
{
class Disk
{
public:
    //Constructor
    Disk() {};      ///< Default constructor
    Disk(Vec3_t const & x0, Vec3_t const & v0, double R0, double M0, double K0, double dt);


    // Methods;
    void StartForce(Vec3_t f = Vec3_t(0.0,0.0,0.0)) {F = f;}
    void DrawDisk(LBM::Lattice & l, double dt);
    void Move(double dt);


    // Data

    Vec3_t x; ///< Position of disk in space
    Vec3_t v; ///< Velocity of the particle
    double R; ///< Radius of the disk
    double M; ///< Mass of the disk
    double K; ///< Disk stiffness
    Vec3_t F; ///< Force over the particle
    Vec3_t xb;///< Position at the previous tme step given by the 


        
};

// Implentation

void CalcForce(Disk * const A,Disk * const B)
{
    double K = 2*(A->K*B->K/(A->K+B->K)); //Equivalent stiffness
    double G = 16;                        //Dissipation constant
    double dist = norm(B->x - A->x);      //Distance between the two spheres
    Vec3_t n = (B->x - A->x)/dist;        //Normal Vector
    double delta = A->R + B->R - dist;    //Overlapping distance
    if (delta>0.0)
    {
        Vec3_t F = K*delta*n;
        Vec3_t vrel = B->v - A->v;
        F -= G*dot(vrel,n)*n;
        A->F -= F;
        B->F += F;
    }
}

Disk::Disk(Vec3_t const & x0, Vec3_t const & v0, double R0, double M0, double K0, double dt)
    : x(x0), v(v0), R(R0), M(M0), K(K0)
{
    xb = x - v*dt;
} 

void Disk::DrawDisk(LBM::Lattice & l, double dt)
{
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		// Set solids and set the magnus velocity
		if (pow(x(0)-double(i),2.0) + pow(x(1)-double(j),2.0) <= pow(R,2.0)) // circle equation
			l.GetCell(i,j)->SetSolid(v(0), v(1));
	}

	 //Calculate de force in the solid
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
				F(0) += (2.0*c->F(op) - alpha*(c->C(op,0)*v(0)+c->C(op,1)*v(1)))/dt*c->C(op,0);
				F(1) += (2.0*c->F(op) - alpha*(c->C(op,0)*v(0)+c->C(op,1)*v(1)))/dt*c->C(op,1);
			}
		}
	}

}

void Disk::Move(double dt)
{
    Vec3_t temp,xa;
    xa    = 2*x - xb + F*(dt*dt/M);
    temp  = xa - x;
    v    = 0.5*(xa - xb)/dt;
    xb   = x;
    x    = xa;
}

}
#endif // MECHSYS_LBM_DEM_H
