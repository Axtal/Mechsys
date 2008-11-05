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
	double errors = 0.0;
	
	/////////////////////////////////////////////////////////////////////////////////////////// 2D - 1 /////
	
	{
		// Set input PSLG (polygon)
		Mesh::Unstructured mu(/*Is3D*/false);

		// Define PSLG polygon sizes
		mu.SetPolySize (/*NPoints*/6, /*NSegments*/5, /*NRegions*/2);

		// Points
		mu.SetPolyPoint (0, 0, 0); // iPoint, X, Y
		mu.SetPolyPoint (1, 1, 0);
		mu.SetPolyPoint (2, 1, 1);
		mu.SetPolyPoint (3, 0, 1);
		mu.SetPolyPoint (4, 0, 0.5);
		mu.SetPolyPoint (5, 1, 0.5);

		// Segments
		mu.SetPolySegment (0, 0, 1, -10); // iSegment, iPointLeft, iPointRight, Tag
		mu.SetPolySegment (1, 1, 2, -20);
		mu.SetPolySegment (2, 2, 3, -30);
		mu.SetPolySegment (3, 3, 0, -40);
		mu.SetPolySegment (4, 4, 5);

		// Regions
		mu.SetPolyRegion (0, -123, 0.01,  0.5, 0.25); // iRegion, Tag, MaxArea, X, Y
		mu.SetPolyRegion (1, -321,  0.1,  0.5, 0.75);

		// Generate
		std::clock_t start = std::clock(); // Initial time
		size_t ne = mu.Generate ();
		std::clock_t total = std::clock() - start; // Time elapsed
		std::cout << "2D("<<ne<<" elements): Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

		cout << mu << endl;

		// Output
		mu.WriteVTU ("tmesh03_2D_1.vtu");
		cout << "[1;34mFile <tmesh03_2D_1.vtu> created[0m" << endl << endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////// 2D - 2 /////

	{
		// Set input PSLG (polygon)
		Mesh::Unstructured mu(/*Is3D*/false);

		// Define PSLG polygon sizes
		mu.SetPolySize (/*NPoints*/29, /*NSegments*/29, /*NRegions*/0, /*NHoles*/1);

		// Points
		mu.SetPolyPoint ( 0, 0.200000, -0.776400); // iPoint, X, Y
		mu.SetPolyPoint ( 1, 0.220000, -0.773200);
		mu.SetPolyPoint ( 2, 0.245600, -0.756400);
		mu.SetPolyPoint ( 3, 0.277600, -0.702000);
		mu.SetPolyPoint ( 4, 0.488800, -0.207600);
		mu.SetPolyPoint ( 5, 0.504800, -0.207600);
		mu.SetPolyPoint ( 6, 0.740800, -0.739600);
		mu.SetPolyPoint ( 7, 0.756000, -0.761200);
		mu.SetPolyPoint ( 8, 0.774400, -0.772400);
		mu.SetPolyPoint ( 9, 0.800000, -0.776400);
		mu.SetPolyPoint (10, 0.800000, -0.792400);
		mu.SetPolyPoint (11, 0.579200, -0.792400);
		mu.SetPolyPoint (12, 0.579200, -0.776400);
		mu.SetPolyPoint (13, 0.621600, -0.771600);
		mu.SetPolyPoint (14, 0.633600, -0.762800);
		mu.SetPolyPoint (15, 0.639200, -0.744400);
		mu.SetPolyPoint (16, 0.620800, -0.684400);
		mu.SetPolyPoint (17, 0.587200, -0.604400);
		mu.SetPolyPoint (18, 0.360800, -0.604400);
		mu.SetPolyPoint (19, 0.319200, -0.706800);
		mu.SetPolyPoint (20, 0.312000, -0.739600);
		mu.SetPolyPoint (21, 0.318400, -0.761200);
		mu.SetPolyPoint (22, 0.334400, -0.771600);
		mu.SetPolyPoint (23, 0.371200, -0.776400);
		mu.SetPolyPoint (24, 0.371200, -0.792400);
		mu.SetPolyPoint (25, 0.374400, -0.570000);
		mu.SetPolyPoint (26, 0.574400, -0.570000);
		mu.SetPolyPoint (27, 0.473600, -0.330800);
		mu.SetPolyPoint (28, 0.200000, -0.792400);

		// Segments
		mu.SetPolySegment ( 0,28, 0, 0); // iSegment, iPointLeft, iPointRight, Tag
		mu.SetPolySegment ( 1, 0, 1, 0);
		mu.SetPolySegment ( 2, 1, 2, 0);
		mu.SetPolySegment ( 3, 2, 3, 0);
		mu.SetPolySegment ( 4, 3, 4, 0);
		mu.SetPolySegment ( 5, 4, 5, 0);
		mu.SetPolySegment ( 6, 5, 6, 0);
		mu.SetPolySegment ( 7, 6, 7, 0);
		mu.SetPolySegment ( 8, 7, 8, 0);
		mu.SetPolySegment ( 9, 8, 9, 0);
		mu.SetPolySegment (10, 9,10, 0);
		mu.SetPolySegment (11,10,11, 0);
		mu.SetPolySegment (12,11,12, 0);
		mu.SetPolySegment (13,12,13, 0);
		mu.SetPolySegment (14,13,14, 0);
		mu.SetPolySegment (15,14,15, 0);
		mu.SetPolySegment (16,15,16, 0);
		mu.SetPolySegment (17,16,17, 0);
		mu.SetPolySegment (18,17,18, 0);
		mu.SetPolySegment (19,18,19, 0);
		mu.SetPolySegment (20,19,20, 0);
		mu.SetPolySegment (21,20,21, 0);
		mu.SetPolySegment (22,21,22, 0);
		mu.SetPolySegment (23,22,23, 0);
		mu.SetPolySegment (24,23,24, 0);
		mu.SetPolySegment (25,24,28, 0);
		mu.SetPolySegment (26,25,26, 0);
		mu.SetPolySegment (27,26,27, 0);
		mu.SetPolySegment (28,27,25, 0);

		// Holes
		mu.SetPolyHole (0, 0.47, -0.5); // iHole, X, Y

		// Generate
		std::clock_t start = std::clock(); // Initial time
		size_t ne = mu.Generate ();
		std::clock_t total = std::clock() - start; // Time elapsed
		std::cout << "2D("<<ne<<" elements): Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

		// Output
		mu.WriteVTU ("tmesh03_2D_2.vtu");
		cout << "[1;34mFile <tmesh03_2D_2.vtu> created[0m" << endl << endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////// 3D /////
	
	{
	}

	// Return error flag
	if (fabs(errors)>1.0e-13) return 1;
	else                      return 0;
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
