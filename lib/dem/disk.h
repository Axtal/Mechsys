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

#ifndef DEM_DISK_H
#define DEM_DISK_H

#include <iostream>
#include <math.h>
#include <string>
#include <vector>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

typedef blitz::TinyVector<double,2> Vec2_t;

class Disk
{
public:
	// Constructors
	Disk () {}
	Disk (double rho0, double r0, Vec2_t x0, Vec2_t v0, double dt0);
	
	// Methods
	void Move();

	// Access methods
	double & rho () { return _rho; } ///< Density
	double & r   () { return _r;   } ///< Radious
	double & dt  () { return _dt;  } ///< time step
	double & V   () { return _V;   } ///< Volume
	double & m   () { return _m;   } ///< mass
	Vec2_t & x   () { return _x;   } ///< _Position
	Vec2_t & v   () { return _v;   } ///< Velocity
	Vec2_t & xp  () { return _xp;  } ///< previous position
	Vec2_t & F   () { return _F;   } ///< Force 

protected:
	//Variables
	double _rho; ///< Density
	double _r;   ///< Radious
	double _dt;  ///< time step
	double _V;   ///< Volume
	double _m;   ///< mass
	Vec2_t _x;   ///< Position
	Vec2_t _v;   ///< Velocity
	Vec2_t _xp;  ///< previous position
	Vec2_t _F;   ///< Force
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Disk::Disk(double rho0, double r0, Vec2_t x0, Vec2_t v0, double dt0)
{
	_rho = rho0;				
	_r   = r0;		
	_x   = x0;		
	_v   = v0;		
	_dt  = dt0;		
	_xp  = _x - _v*_dt;	
	_V   = M_PI*_r*_r;	
	_m   = _V*_rho;		
	_F   = 0.0,0.0;		
}

inline void Disk::Move()
{
	Vec2_t tmp = _x;
	_x         = 2*_x-_xp + (_F/_m)*_dt*_dt;
	_v         = (_x - _xp)/(2*_dt);
	_xp        = tmp;
}

#endif // DEM_DISK_H
