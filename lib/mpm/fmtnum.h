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

/* FmtNum - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_FMTNUM_H
#define MPM_FMTNUM_H

// STL
#include <iostream>
#include <iomanip>

namespace MPM {

/** Formatting Numbers via STL streams. */
struct FmtNum
{
	bool BoolAlpha;  ///< Format output as a boolean?
	bool Integer;    ///< Format output as an integer?
	bool Scientific; ///< Format output as a scientific number?
	int  Width;      ///< Width of the output
	int  Precision;  ///< Precision for floating point numbers
};

//                 bool  integ  scien   w   p
FmtNum _a     = {  true, false, false,  0,  0 }; ///< Boolean
FmtNum _3     = { false,  true, false,  3,  0 }; ///< Integer
FmtNum _4     = { false,  true, false,  4,  0 }; ///< Integer
FmtNum _6     = { false,  true, false,  6,  0 }; ///< Integer
FmtNum _8     = { false,  true, false,  8,  0 }; ///< Integer
FmtNum _3s    = { false, false,  true,  0,  3 }; ///< Scientific
FmtNum _8s    = { false, false,  true,  0,  8 }; ///< Scientific
FmtNum _4_2   = { false, false, false,  4,  2 }; ///< General
FmtNum _6_3   = { false, false, false,  6,  3 }; ///< General
FmtNum _12_6  = { false, false, false, 12,  6 }; ///< General
FmtNum _20_15 = { false, false, false, 20, 15 }; ///< General

/** Format the output. */
std::ostream & operator<< (std::ostream & os, FmtNum const & FN)
{
	     if (FN.BoolAlpha)  { os<<" "<<std::setw(6)<<std::boolalpha; return os; }
	else if (FN.Integer)    { os<<" "<<std::setw(FN.Width); return os; }
	else if (FN.Scientific) { os<<" "<<std::setw(FN.Precision+9)<<std::scientific<<std::setprecision(FN.Precision); return os; } // add 9
	else                    { os<<" "<<std::setw(FN.Width+1)<<std::fixed<<std::setprecision(FN.Precision); return os; } // add 1 (for sign)
}

}; // namespace MPM

#endif // MPM_FMTNUM_H
