/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

#ifndef MECHSYS_NUMSTREAMS_H
#define MECHSYS_NUMSTREAMS_H

// STL
#include <iostream>
#include <iomanip>

namespace Util
{

/** Number format via STL streams. */
struct NumStream
{
	bool BoolAlpha;  ///< Format output as a boolean?
	bool Integer;    ///< Format output as an integer?
	bool Scientific; ///< Format output as a scientific number?
	int  Width;      ///< Width of the output
	int  Precision;  ///< Precision for floating point numbers
};

//                    bool  integ  scien   w   p
NumStream _a     = {  true, false, false,  0,  0 }; ///< Boolean
NumStream _3     = { false,  true, false,  3,  0 }; ///< Integer
NumStream _4     = { false,  true, false,  4,  0 }; ///< Integer
NumStream _6     = { false,  true, false,  6,  0 }; ///< Integer
NumStream _8     = { false,  true, false,  8,  0 }; ///< Integer
NumStream _3s    = { false, false,  true,  0,  3 }; ///< Scientific
NumStream _8s    = { false, false,  true,  0,  8 }; ///< Scientific
NumStream _6_2   = { false, false, false,  6,  3 }; ///< General
NumStream _6_3   = { false, false, false,  6,  3 }; ///< General
NumStream _6_4   = { false, false, false,  6,  4 }; ///< General
NumStream _6_6   = { false, false, false,  6,  6 }; ///< General
NumStream _8_2   = { false, false, false,  8,  2 }; ///< General
NumStream _8_3   = { false, false, false,  8,  2 }; ///< General
NumStream _8_4   = { false, false, false,  8,  4 }; ///< General
NumStream _12_6  = { false, false, false, 12,  6 }; ///< General
NumStream _20_15 = { false, false, false, 20, 15 }; ///< General

/** Format the output. */
std::ostream & operator<< (std::ostream & os, NumStream const & NS)
{
	     if (NS.BoolAlpha)  { os<<" "<<std::setw(6)<<std::boolalpha; return os; }
	else if (NS.Integer)    { os<<" "<<std::setw(NS.Width); return os; }
	else if (NS.Scientific) { os<<" "<<std::setw(NS.Precision+9)<<std::scientific<<std::setprecision(NS.Precision); return os; } // add 9
	else                    { os<<" "<<std::setw(NS.Width+1)<<std::fixed<<std::setprecision(NS.Precision); return os; } // add 1 (for sign)
}

}; // namespace Util

#endif // MECHSYS_NUMSTREAMS_H
