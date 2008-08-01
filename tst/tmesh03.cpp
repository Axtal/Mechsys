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
#include <ctime>  // for std::clock()

// MechSys
#include "mesh/unstructured.h"
#include "util/array.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "linalg/matrix.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;

int main(int argc, char **argv) try
{
	
	/////////////////////////////////////////////////////////////////////////////////////////// 2D /////
	
	{
		// Set geometry
		Mesh::Unstructured mu;
		mu.SetNVerts (4);
		mu.SetVert   (0, 0,0);
		mu.SetVert   (1, 1,0);
		mu.SetVert   (2, 1,1);
		mu.SetVert   (3, 0,1);

		cout << mu << endl;

		// Generate
		std::clock_t start = std::clock(); // Initial time
		size_t ne = mu.Generate ();
		std::clock_t total = std::clock() - start; // Time elapsed
		std::cout << "2D("<<ne<<" elements): Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

		// Output
		mu.WriteVTU ("tmesh03_2D.vtu");
		cout << "[1;34mFile <tmesh03_2D.vtu> created[0m" << endl << endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////// 3D /////
	
	{
	}

	return 0;
}
catch (Exception * e) 
{
	e->Cout();
	if (e->IsFatal()) {delete e; exit(1);}
	delete e;
}
catch (char const * m)
{
	std::cout << "Fatal: " << m << std::endl;
	exit (1);
}
catch (...)
{
	std::cout << "Some exception (...) ocurred\n";
} 
