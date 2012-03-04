/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso                                *
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

/* CHull2D - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_CHULL2D_H
#define MPM_CHULL2D_H

// STL
#include<algorithm> // for std::sort

// Local
#include <mechsys/mpm/defs.h>
#include <mechsys/util/array.h>

namespace MPM {

/** Lexicographically comparision between points. */
inline bool CompareVec2D (Vector3D const & A, Vector3D const & B)
{
	return (A(0)<B(0) || (A(0)==B(0) && A(1)<B(1)));
}

/** 2D cross product.
 * Return a positive value, if OAB makes a counter-clockwise turn,
 * negative for clockwise turn, and zero if the points are collinear. */
inline double Cross2D (Vector3D const & O, Vector3D const & A, Vector3D const & B)
{
	return (A(0)-O(0)) * (B(1)-O(1)) - (A(1)-O(1)) * (B(0)-O(0));
}

/** Convex hull using Andrew's monotone chain 2D algorithm.
 * Returns a list of points on the convex hull in counter-clockwise order. */
inline void CHull2D (Array<Vector3D> & P, Array<size_t> & H)
{
	// Input
	int n = P.Size();
	int k = 0;
	Array<Vector3D> h(2*n); // points on the hull
	Array<size_t>   j(2*n); // indexes of the points in the hull

	// Sort points lexicographically
	std::sort (P.GetPtr(), P.GetPtr()+P.Size(), CompareVec2D);

	// Build lower hull
	for (int i=0; i<n; i++)
	{
		while (k>=2 && Cross2D (h[k-2], h[k-1], P[i])<=0) k--;
		h[k++] = P[i];
		j[k-1] = i;
	}

	// Build upper hull
	for (int i=n-2, t=k+1; i>=0; i--)
	{
		while (k >= t && Cross2D (h[k-2], h[k-1], P[i])<=0) k--;
		h[k++] = P[i];
		j[k-1] = i;
	}

	// Results
	if (k>1)
	{
		H.Resize(k-1);
		for (int i=0; i<k-1; ++i) H[i] = j[i];
	}
}

}; // namespace MPM

#endif // MPM_CHULL2D_H
