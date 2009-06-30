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
#include <math.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/quaternion.h"

typedef blitz::TinyVector<double,3> Vec3_t;

class Particle
{
public:
	// Constructor
	Particle(double         rho0,    ///< Density
	         Vec3_t       & r0,      ///< Initial position
	         Vec3_t       & v0,      ///< Initial velocity
	         Vec3_t       & omega0,  ///< Initial angular velocity
	         Quaternion_t & q0);     ///< Initial quaternion giving the orientation
	
	// Destructor
	virtual ~Particle() {}

	// Set methods
	void Start              (double   dt);  ///< Start the particle for the Verlet algorithm
	void Translation        (double   dt);  ///< Translates the particle after the forces are calculated
	void Translation        (Vec3_t   t);   ///< Translates the particle a vector t
	void Rotation           (double   dt);  ///< Rotates the particle after the torques are calculated
	void QuaternionRotation (double   dt);  ///< Rotates the particle by a quaternion
	void StartForce         ();             ///< Starts the force at a given value
	void StartForce         (Vec3_t & F);   ///< Starts the force at a given value

	// Special
	virtual void Special() =0; // pure virtual function => derived class MUST derive
	virtual void Another() {}  // not pure => derive clas MAY drive

protected:
	double       _m;    ///< Mass of the particle
	double       _V;    ///< Volume of the particle
	double       _rho;  ///< Density of the particle
	Vec3_t       _r;    ///< Position of the particle
	Vec3_t       _rp;   ///< Previous position of the particle
	Vec3_t       _v;    ///< Velocity of the particle
	Vec3_t       _f;    ///< Force of the particle
	Vec3_t       _T;    ///< Torque of the particle
	Vec3_t       _ome;  ///< Angular velocity of the particle
	Vec3_t       _omep; ///< Previous angular velocity of the particle
	Quaternion_t _q;    ///< Quaternion of the particle
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void Particle::Start(double dt)
{
}

inline void Particle::Translation(double dt)
{
}

inline void Particle::Translation(Vec3_t t)
{
}

inline void Particle::Rotation(double dt)
{
}

inline void Particle::QuaternionRotation(double dt)
{
}

inline void Particle::StartForce()
{
}

inline void Particle::StartForce(Vec3_t & F)
{
}


#endif // DEM_PARTICLE_H
