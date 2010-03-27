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

// STL
#include <iostream>
#include <fstream>

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/unstructured.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	Array<double> a(4);
	a = 1,2,3,4;

	cout << a << endl;

    //std::ofstream ofs("tserialize.out");

    // save data to archive
	//boost::archive::text_oarchive oa(ofs);
	//oa << a;

	return 0;
}
MECHSYS_CATCH
