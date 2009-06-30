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

#ifndef DEM_SPHERE_H
#define DEM_SPHERE_H

// Std lib
#include <math.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/particle.h"
#include "dem/quaternion.h"

typedef blitz::TinyVector<double,3> Vec3_t;

class Sphere : public Particle
{
public:
	// Constructor
	Sphere(double         rho0,    ///< Density
	       Vec3_t       & r0,      ///< Initial position
	       Vec3_t       & v0,      ///< Initial velocity
	       Vec3_t       & omega0,  ///< Initial angular velocity
	       Quaternion_t & q0,      ///< Initial quaternion giving the orientation
	       double         R0);     ///< Radius
	
	void Special() {}

protected:
	double _R; ///< Radius of the particle
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructor
inline Sphere::Sphere(double rho0, Vec3_t & r0, Vec3_t & v0, Vec3_t & omega0, Quaternion_t & q0, double R0)
	: Particle(rho0,r0,v0,omega0,q0),
	  _R (R0)
{
}

#endif // DEM_SPHERE_H
