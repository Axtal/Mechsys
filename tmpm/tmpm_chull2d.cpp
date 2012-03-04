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

/* tchull - Copyright (C) 2007 Dorival de Moraes Pedroso */

// STL
#include <iostream>

// MechSys
#include <mechsys/util/fatal.h>

// Local
#include <mechsys/mpm/defs.h>
#include <mechsys/mpm/fmtnum.h>
#include <mechsys/mpm/chull2d.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	Array<MPM::Vector3D> Ps(5);
	Ps[0] = 1,1;
	Ps[1] = 1,0;
	Ps[2] = 0,1;
	Ps[3] = 0,0;
	Ps[4] = 0.5,0.5;

	Array<size_t> Hs;
    MPM::CHull2D (Ps, Hs);
	cout << "Ps = \n" << Ps << endl;
	cout << "Hs = \n" << Hs << endl;

	return 0;
}
MECHSYS_CATCH
