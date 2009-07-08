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

#ifndef DEM_FACE3D_H
#define DEM_FACE3D_H

// Std lib
#include <math.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/edge3d.h"
#include "util/array.h"

class Face3D
{
public:
	// Constructor
	Face3D(void) {};          ///< Default Constructor
	Face3D(const Vec3_t * a,  ///< Vector array
               const size_t N);   ///< Numer of sides

	// Access Methods
	Edge3D * Edge (size_t i) {return _sides[i];}          ///< Returns pointer to the i-th side
	size_t NumberofSides () {return (int) _sides.Size();} ///< Returns the number of sides



protected:
	Array<Edge3D *> _sides;   ///< Array of edges representing sides of the face.
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Face3D::Face3D (const Vec3_t * a,const size_t N)
{
	_sides.Resize (N);
	for(size_t i=0;i<N;i++) 
	{
		_sides[i] = new Edge3D(a[(i+1)%N],a[i]);
	}
}

#endif //DEM_FACE3D_H
