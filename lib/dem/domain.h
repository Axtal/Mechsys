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
#include "dem/interacton.h"
#include "util/array.h"

class Domain
{
public:
	//Constructor
	Domain() {}

	//Methods
	void GenerateSpheres(size_t n,                                 ///< Number of spheres
	                     double Xmin,                              ///< Left boundary
	                     double Xmax,                              ///< Right boundary
	                     double Ymin,                              ///< Back boundary
	                     double Ymax,                              ///< Front boundary
	                     double Zmin,                              ///< Bottom boundary
	                     double Zmax,                              ///< Top boundary
	                     double rho,                               ///< Density of the material
	                     double Rmin);                             ///< Minimun radius in units of the maximun radius
	
	void GenerateBox(double Xmin,                                  ///< Left boundary
	                 double Xmax,                                  ///< Right boundary
	                 double Ymin,                                  ///< Back boundary
	                 double Ymax,                                  ///< Front boundary
	                 double Zmin,                                  ///< Bottom boundary
	                 double Zmax,                                  ///< Top boundary
	                 double Thickness);                            ///< Thickness of the wall, cannot be zero

	void AddTetra(const Vec3_t & r,double R,double l,double rho0); ///< Add a tetrahedron at position r with spheroradius R, side of length l and density rho0
	void AddRice(const Vec3_t & r,double R,double l,double rho0);  ///< Add a rice at position r with spheroradius R, side of length l and density rho0
	void AddCube(const Vec3_t & r,double R,double l,double rho0);  ///< Add a cube at position r with spheroradius R, side of length l and density rho0

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
	size_t nx = Lx/(2*Rmax),ny = Ly/(2*Rmax);
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

inline void Domain::GenerateBox(double Xmin,double Xmax,double Ymin,double Ymax,double Zmin,double Zmax,double Thickness)
{
	
}

inline void Domain::AddTetra(const Vec3_t & r,double R,double l,double rho0)
{
	Array<Vec3_t *> V;
	V.Push(new Vec3_t(l/sqrt(8),l/sqrt(8),l/sqrt(8)));
	V.Push(new Vec3_t(-l/sqrt(8),-l/sqrt(8),l/sqrt(8)));
	V.Push(new Vec3_t(-l/sqrt(8),l/sqrt(8),-l/sqrt(8)));
	V.Push(new Vec3_t(l/sqrt(8),-l/sqrt(8),-l/sqrt(8)));
	Array<Array <int> > E;
	E.Resize(6);
	size_t n = 0;
	for (size_t i = 0;i < 3;i++)
	{
		for (size_t j = i+1;j < 4;j++)
		{
			E[n].Resize(2);
			E[n][0] = i;
			E[n][1] = j;
			n++;
		}
	}
	Array<Array <int> > F;
	F.Resize(4);
	F[0].Push(0);
	F[0].Push(1);
	F[0].Push(2);

	F[1].Push(0);
	F[1].Push(1);
	F[1].Push(3);

	F[2].Push(0);
	F[2].Push(2);
	F[2].Push(3);

	F[3].Push(1);
	F[3].Push(2);
	F[3].Push(3);
	double angle = (1.*rand())/RAND_MAX*2*M_PI;
	Vec3_t axis;
	axis(0) = (1.*rand())/RAND_MAX;
	axis(1) = (1.*rand())/RAND_MAX;
	axis(2) = (1.*rand())/RAND_MAX;
	Quaternion_t q;
	NormalizeRotation(angle,axis,q);

	for (size_t i = 0;i < V.Size();i++)
	{
		Vec3_t t;
		Rotation(*V[i],q,t);
		*V[i]=t+r;
	}
	Vec3_t Zero(0,0,0);
	_Particles.Push(new Particle(V,E,F,R,rho0,Zero,Zero));
}

inline void Domain::AddRice(const Vec3_t & r,double R,double l,double rho0)
{
	Array<Vec3_t *> V;
	V.Push(new Vec3_t(0,0,l/2));
	V.Push(new Vec3_t(0,0,-l/2));
	Array<Array <int> > E;
	E.Resize(1);
	E[0].Resize(2);
	E[0][0] = 0;
	E[0][1] = 1;
	Array<Array <int> > F;
	F.Resize(0);
	double angle = (1.*rand())/RAND_MAX*2*M_PI;
	Vec3_t axis;
	axis(0) = (1.*rand())/RAND_MAX;
	axis(1) = (1.*rand())/RAND_MAX;
	axis(2) = (1.*rand())/RAND_MAX;
	Quaternion_t q;
	NormalizeRotation(angle,axis,q);
	//q=1,0,0,0;
	std::cout<<q<<std::endl;
	for (size_t i = 0;i < V.Size();i++)
	{
		Vec3_t t;
		Rotation(*V[i],q,t);
		*V[i]=t+r;
	}
	Vec3_t Zero(0,0,0);
	_Particles.Push(new Particle(V,E,F,R,rho0,Zero,Zero));
}

inline void Domain::AddCube(const Vec3_t & r,double R,double l,double rho0)
{
	Array<Vec3_t *> V;
	for (size_t i = 0;i < 4;i++)
	{
		V.Push(new Vec3_t(l*cos((2*i+1)*M_PI/4)/sqrt(2),l*sin((2*i+1)*M_PI/4)/sqrt(2),l/2));
		V.Push(new Vec3_t(l*cos((2*i+1)*M_PI/4)/sqrt(2),l*sin((2*i+1)*M_PI/4)/sqrt(2),-l/2));
	}
	Array<Array <int> > E;
	E.Resize(12);
	for (size_t i = 0;i < 4;i++)
	{
		E[i].Push(2*i);
		E[i].Push(2*i+1);
		E[i+4].Push((2*i)%8);
		E[i+4].Push((2*i+2)%8);
		E[i+8].Push((2*i+1)%8);
		E[i+8].Push((2*i+3)%8);
	}
	Array<Array <int> > F;
	F.Resize(6);
	for (size_t i = 0;i < 6;i++) \
	{
		F[i].Resize(4);
	}
	F[0][0] = 0; F[1][0] = 1; F[2][0] = 0; F[3][0] = 6; F[4][0] = 2; F[5][0] = 0;
	F[0][1] = 2; F[1][1] = 3; F[2][1] = 2; F[3][1] = 4; F[4][1] = 4; F[5][1] = 6;
	F[0][2] = 4; F[1][2] = 5; F[2][2] = 3; F[3][2] = 5; F[4][2] = 5; F[5][2] = 7;
	F[0][3] = 6; F[1][3] = 7; F[2][3] = 1; F[3][3] = 7; F[4][3] = 3; F[5][3] = 1;
	double angle = (1.*rand())/RAND_MAX*2*M_PI;
	Vec3_t axis;
	axis(0) = (1.*rand())/RAND_MAX;
	axis(1) = (1.*rand())/RAND_MAX;
	axis(2) = (1.*rand())/RAND_MAX;
	Quaternion_t q;
	NormalizeRotation(angle,axis,q);

	for (size_t i = 0;i < V.Size();i++)
	{
		Vec3_t t;
		Rotation(*V[i],q,t);
		*V[i]=t+r;
	}
	Vec3_t Zero(0,0,0);
	_Particles.Push(new Particle(V,E,F,R,rho0,Zero,Zero));
}


#endif //DEM_DOMAIN_H
