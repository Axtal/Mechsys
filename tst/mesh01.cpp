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

// MechSys
#include "util/exception.h"
#include "mesh/unstructured.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	Mesh::Unstructured mesh(/*Is3D*/false);
	mesh.SetPolySize    (/*NPoints*/4, /*NSegments*/4, /*NRegions*/1);
	mesh.SetPolyPoint   (0, 0.0, 0.0);
	mesh.SetPolyPoint   (1, 1.0, 0.0);
	mesh.SetPolyPoint   (2, 1.0, 1.0);
	mesh.SetPolyPoint   (3, 0.0, 1.0);
	mesh.SetPolySegment (0, 0, 1);
	mesh.SetPolySegment (1, 1, 2);
	mesh.SetPolySegment (2, 2, 3);
	mesh.SetPolySegment (3, 3, 0);
	mesh.SetPolyRegion  (0, /*Tag*/-1, /*MaxArea*/-1, /*Xc*/0.5, /*Yc*/0.5);
	mesh.Generate       (/*WithInfo*/true);
	mesh.WriteVTU       ("tmesh01.vtu");

	return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
