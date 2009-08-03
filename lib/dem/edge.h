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

#ifndef DEM_EDGE3D_H
#define DEM_EDGE3D_H

// Std lib
#include <math.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/quaternion.h"

class Edge
{
public:
	// Constructor
	Edge(void) {};          ///< Default Constructor
	Edge(const Vec3_t & a,  ///< Initial vector
	       const Vec3_t & b); ///< Final vector

	
	// Access Methods
	Vec3_t & ri () {return _ri;} ///< Initial vector
	Vec3_t & rf () {return _rf;} ///< Final vector
	Vec3_t & dr () {return _dr;} ///< Difference vector

	// Methods
	void Rotate(const Quaternion_t & q, ///< Quaternion representing the rotation
                    const Vec3_t & v);      ///< Position of the axis of rotation



protected:
	double _l;  ///< Length of the Edge
	Vec3_t _ri; ///< Initial position
	Vec3_t _rf; ///< Final position
	Vec3_t _dr; ///< Difference Vector
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Edge::Edge (const Vec3_t & a,const Vec3_t & b)
{
	_ri = a;
	_rf = b;
	_dr = _rf-_ri;
	_l  = norm(_dr);
}



inline void Edge::Rotate (const Quaternion_t & q,const Vec3_t & v)
{
	Vec3_t t1,t2;
	t1 = _ri - v;
	t2 = _rf - v;
	Rotation(t1,q,_ri);
	Rotation(t2,q,_rf);
	_ri = _ri + v;
	_rf = _rf + v;
	_dr = _rf - _ri;
}

#endif //DEM_EDGE3D_H
