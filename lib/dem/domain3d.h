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
#include <map>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "util/array.h"
#include "dem/particle.h"
#include "dem/sphere.h"
#include "dem/interacton.h"
#include "dem/quaternion.h"

class Domain3D
{
	Domain3D();
    ~Domain3D();
protected:
	Array<Particle *> _particles; ///< Array with all particles
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Domain3D::Domain3D()
{
	Vec3_t       r0,v0,omega0;
	Quaternion_t q0;
	double       rho0 = 1.0;
	double       R0   = 1.0;
	r0     = 0.0, 0.0, 0.0;
	v0     = 1.0, 1.0, 1.0;
	omega0 = 1.0, 1.0, 1.0;
	q0     = 1.0, 1.0, 1.0, 1.0;

	Sphere s(rho0, r0, v0, omega0, q0, R0);

	//for (size_t i=0; i<10; ++i)
		//_particles.Push(new Sphere(rho0, r0, v0, omega0, q0, R0));
}

inline Domain3D::~Domain3D()
{
	for (size_t i=0; i<_particles.Size(); ++i)
		if (_particles[i]!=NULL) delete _particles[i];
}

#endif // DEM_DOMAIN_H
