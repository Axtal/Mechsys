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

#ifndef DEM_PARTICLE_H
#define DEM_PARTICLE_H

// Std lib
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/edge3d.h"
#include "util/array.h"

class Particle
{
public:
	// Constructors
	// Sphere constructor, the simplest one
	Particle(const double Radius,     ///< Radius for the sphere
                 const double rho0,       ///< Density of mass
	         const Vec3_t & r0,       ///< Initial position
	         const Vec3_t & v0,       ///< Initial velocity
	         const Vec3_t & ome0)     ///< Initial angular velocity

	// Plane constructor, for walls and the like
	Particle(const double thickness,  ///< Thickness of the plane, equivalent to the spheroradius
                 const double rho0,       ///< Density of mass
	         const Vec3_t & r0,       ///< Initial position
	         const Vec3_t & v0,       ///< Initial velocity
	         const Vec3_t & ome0,     ///< Initial angular velocity
	         const Quaternion_t q0)   ///< Initial orientation
	         


protected:

	bool            _IsLarge; ///< Flag to see if it is a large particle
	double          _m;       ///< Mass of the particle
	double          _R;       ///< Spheroradius of the particle
	Vec3_t          _r;       ///< Position of the particle
	Vec3_t          _rb;      ///< Past position for the Verlet algorithm
	Vec3_t          _v;       ///< Velocity of the particle
	Vec3_t          _f;       ///< Force over the particle
	Vec3_t          _T;       ///< Torque over the particle
	Vec3_t          _w;       ///< Angular velocity
	Vec3_t          _wb;      ///< Past angular velocity for the leap frog algorithm
	Vec3_t          _I;       ///< Vector containing the principal components of the inertia tensor
	Quaternion_t    _Q;       ///< The quaternion representing the rotation
	Array<Vec3_t *> _vertex;  ///< The set of vertices defining the geometry of the particle
	Array<Edge3D *> _edges;   ///< The set of edges defining the geometry of the particle
	Array<Face3D *> _faces;   ///< The set of faces defining the geometry of the particle
	
	
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////



#endif // DEM_PARTICLE_H
