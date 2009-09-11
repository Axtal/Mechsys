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

inline bool IsNan(double Val)
{
    return (std::isnan(Val) || ((Val==Val)==false)); // NaN is the only value, for which the expression Val==Val is always false
}

/*
inline double Sgn        (double Val)              { return (Val>=0.0 ? +1.0 : -1.0);                                } ///< Sgn function where Sgn(0)=+1
inline double Signal     (double Val, double Zero) { return (fabs(Val)<=Zero ? 0.0 : Sgn(Val));                      } ///< Sgn function where Sgn(0)=0
inline double Sign       (double a, double b)      { return (b>=0.0 ? fabs(a) : -fabs(a));                           } ///< Composite Sgn function. Returns |a| or -|a| according to the sign of b
inline double Acos       (double Val)              { return (Val>=1.0 ?  0.0 : (Val<=-1.0 ? Pi() : acos(Val)) );     } ///< Safe acos function
inline bool   Str2Bool   (String const & Str)      { if (Str==String("true") || Str==String("TRUE") || Str==String("True")) return true; else return false; } ///< Converts "TRUE", "True", or "true" to bool
inline bool   IsNanOrInf (double Val)              { int r=std::fpclassify(Val); if (r==FP_NAN) return true; if (r==FP_INFINITE) return true; return false; } ///< Check whether a number is NaN of Inf
*/

template <typename Val_T> inline Val_T Min (Val_T const & a, Val_T const & b) { return (a<b ? a : b); } ///< Minimum between a and b
template <typename Val_T> inline Val_T Max (Val_T const & a, Val_T const & b) { return (a>b ? a : b); } ///< Maximum between a and b

/** Swap two values. */
template <typename Val_T>
inline void Swap (Val_T & a, Val_T & b)
{
    Val_T tmp = a;
    a = b;
    b = tmp;
}

/** Sort two values on an ascending order. */
template <typename Val_T>
inline void Sort (Val_T & a, Val_T & b)
{
    if (b<a) Util::Swap (a,b);
}

/** Sort three values on an ascending order. */
template <typename Val_T>
inline void Sort (Val_T & a, Val_T & b, Val_T & c)
{
    if (b<a) Util::Swap (a,b);
    if (c<b) Util::Swap (b,c);
    if (b<a) Util::Swap (a,b);
}

/** Find best square for given rows and columns. */
/*
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
*/

}; // namespace Util

#endif // MECHSYS_SORT_H
