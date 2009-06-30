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

#ifndef MECHSYS_DEM3D_SPHERE_H
#define MECHSYS_DEM3D_SPHERE_H

#include <math.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

#include "dem3D/vec3.h"

using namespace DEM3D;

class Sphere
{
public:
	Sphere(	double rho0,	///< Density
	       	double R0,	///< Radious
		Vec3_t &r0,	///< Initial position
		Vec3_t &v0,	///< Initial velocity
		Vec3_t &omega0,	///< Initial angular velocity
		Quaternion &q0	///< Initial quaternion giving the orientation
	);
	
	void start(double dt);               ///< It start the particle for the Verlet algorithm
	void translation(double dt);         ///< Translates the particle after the forces are calculated
	void translation(Vec3_t t);          ///< Translates the particle a vector t
	void rotation(double dt);            ///< Rotates the particle after the torques are calculated
	void rotation_quaternion(double dt); ///< Rotates the particle by a quaternion
	void startforce(Vec3_t & F);         ///< Starts the force at a given value

protected:
	double _m;    ///< Mass of the particle
	double _V;    ///< Volume of the particle
	double _rho;  ///< Density of the particle
	double _R;    ///< Radious of the particle
	Vec3_t _r;    ///< Position of the particle
	Vec3_t _rp;   ///< Previous position of the particle
	Vec3_t _v;    ///< Velocity of the particle
	Vec3_t _f;    ///< Force of the particle
	Vec3_t _T;    ///< Torque of the particle
	Vec3_t _ome;  ///< Angular velocity of the particle
	Vec3_t _omep; ///< Previous angular velocity of the particle
	Quaternion q; ///< Quaternion of the particle

};

#endif

