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

#ifndef DEM_DOMAIN_H
#define DEM_DOMAIN_H

// Std lib
#include <math.h>
#include <stdlib.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/particle.h"
#include "util/array.h"

class Domain
{
public:
	//Constructor
	Domain() {}

	//Methods
	void GenerateSpheres(size_t n,                  ///< Number of spheres
	                     double Xmin,               ///< Left boundary
	                     double Xmax,               ///< Right boundary
	                     double Ymin,               ///< Back boundary
	                     double Ymax,               ///< Front boundary
	                     double Zmin,               ///< Bottom boundary
	                     double Zmax,               ///< Top boundary
	                     double rho,                ///< Density of the material
	                     double Rmin);              ///< Minimun radius in units of the maximun radius
	
	void GenerateBox(double Xmin,                   ///< Left boundary
	                 double Xmax,                   ///< Right boundary
	                 double Ymin,                   ///< Back boundary
	                 double Ymax,                   ///< Front boundary
	                 double Zmin,                   ///< Bottom boundary
	                 double Zmax,                   ///< Top boundary
	                 double Thickness);             ///< Thickness of the wall, cannot be zero

	//Access Methods
        size_t NumberParticles ( )         { return   _Particles.Size();} ///< Return the number of particles
        Particle * Particles   ( size_t i) { return _Particles[i];}       ///< Return pointer to i-th particle

protected:
	Array<Particle *> _Particles;                   ///< The array containing all the particles in the Domain.

};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline void Domain::GenerateSpheres(size_t n,double Xmin,double Xmax,double Ymin,double Ymax,double Zmin,double Zmax,double rho,double Rmin)
{
	double Lx = Xmax-Xmin,Ly = Ymax-Ymin,Lz = Zmax-Zmin;
	double Rmax = pow(Lx*Ly*Lz/(8*n),1./3.);
	size_t nx = Lx/(2*Rmax),ny = Ly/(2*Rmax),nz = Lz/(2*Rmax);
	Array <Vec3_t *> Vertex;
	Array <Array <int> > Empty;
	Vertex.Resize(1);
	Empty.Resize(0);
	Vec3_t Zero(0,0,0);
	for (size_t i=0;i<n;i++)
	{
		Vec3_t r(Xmin + Rmax+2*Rmax*(i%nx),Ymin + Rmax+2*Rmax*((i/nx)%ny),Zmin + Rmax+2*Rmax*(i/(nx*ny)));
		double R = Rmax*Rmin + (1.*rand())/RAND_MAX*Rmax*(1-Rmin);
		Vertex[0] = new Vec3_t(0,0,0);
		*Vertex[0] = r;
		Particle *nparticle = new Particle(Vertex,Empty,Empty,R,rho,Zero,Zero);
		_Particles.Push(nparticle);
	}
}




#endif //DEM_DOMAIN_H
