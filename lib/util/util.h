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

#ifndef MECHSYS_SORT_H
#define MECHSYS_SORT_H

// STL
#include <cmath>
#include <cfloat> // for DBL_EPSILON

// MechSys
#include "util/string.h"

namespace Util
{

// Constants
const double ZERO   = sqrt(DBL_EPSILON); ///< Machine epsilon (smaller positive)
const double SQ2    = sqrt(2.0);         ///< \f$ \sqrt{2} \f$
const double SQ3    = sqrt(3.0);         ///< \f$ \sqrt{3} \f$
const double SQ6    = sqrt(6.0);         ///< \f$ \sqrt{6} \f$
const double SQ2BY3 = sqrt(2.0/3.0);     ///< \f$ \sqrt{2/3} \f$
const double PI     = 4.0*atan(1.0);     ///< \f$ \pi \f$

/*
inline double Pi         ()                        { return 3.14159265358979323846264338327950288419716939937510582; } ///< The constant PI
inline double ToRad      (double deg_angle)        { return 0.0174532925199433*deg_angle;                            } ///< Converts degrees to radians
inline double ToDeg      (double rad_angle)        { return 57.2957795130823*rad_angle;                              } ///< Converts radians to degrees
inline double Sgn        (double Val)              { return (Val>=0.0 ? +1.0 : -1.0);                                } ///< Sgn function where Sgn(0)=+1
inline double Signal     (double Val, double Zero) { return (fabs(Val)<=Zero ? 0.0 : Sgn(Val));                      } ///< Sgn function where Sgn(0)=0
inline double Sign       (double a, double b)      { return (b>=0.0 ? fabs(a) : -fabs(a));                           } ///< Composite Sgn function. Returns |a| or -|a| according to the sign of b
inline double Acos       (double Val)              { return (Val>=1.0 ?  0.0 : (Val<=-1.0 ? Pi() : acos(Val)) );     } ///< Safe acos function
inline bool   Str2Bool   (String const & Str)      { if (Str==String("true") || Str==String("TRUE") || Str==String("True")) return true; else return false; } ///< Converts "TRUE", "True", or "true" to bool
inline bool   IsNanOrInf (double Val)              { int r=std::fpclassify(Val); if (r==FP_NAN) return true; if (r==FP_INFINITE) return true; return false; } ///< Check whether a number is NaN of Inf
inline double Sq2        ()                        { return 1.41421356237310;              } ///< Square root of 2.0
inline double Sq3        ()                        { return 1.73205080756888;              } ///< Square root of 3.0
inline double Sq6        ()                        { return 2.44948974278318;              } ///< Square root of 6.0
template <typename Type>
inline Type   Min        (Type a, Type b)          { return (a<b ? a : b); } ///< Minimum between a and b
template <typename Type>
inline Type   Max        (Type a, Type b)          { return (a>b ? a : b); } ///< Maximum between a and b
*/

/** Swap two values. */
inline void Swap(double & x, double & y)
{
    double temp = x;
    x = y;
    y = temp;
}

/** Sort an array according to an ascending order. */
inline void Sort(double A[], int Size)
{
	for (int i=0; i<Size-1; ++i)
	{
		int min_index = i;

		// Find the index of the minimum element
		for (int j=i+1; j<Size; ++j)
			if (A[j] < A[min_index])
				min_index = j;

		// Swap if i-th element not already smallest
		if (min_index > i)
			Util::Swap(A[i], A[min_index]);
	}
}

/** Find best square for given rows and columns. */
inline void FindBestSquare (int Size, int & nRow, int & nCol)
{
	nRow = -1;  // not found
	nCol = -1;  // not found
	for (int x=1; x<=Size; ++x)
	{
		if ((x*x)>=Size)
		{
			if ((x*x)==Size)
			{
				nRow = x;
				nCol = x;
				return;
			}
			else
			{
				for (int y=x; y>=1; --y)
				{
					if ((x*y)==Size)
					{
						nRow = x;
						nCol = y;
						return;
					}
				}
			}
		}
	}
}

}; // namespace Util

#endif // MECHSYS_SORT_H
