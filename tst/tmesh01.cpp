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
#include "mesh/structured.h"
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
		// Constants
		double L     = 1.0;   // length
		double H     = 1.0;   // height
		double hL    = L/2.0; // half-length
		double hH    = H/2.0; // half-height
		size_t ndivx = 30;    // divisions along x
		size_t ndivy = 40;    // divisions along y

		// Coordinates
		Matrix<double> c(2,8);
		c = 0.,  L, L, 0.,    hL,  L, hL, 0.,
			0., 0., H,  H,    0., hH,  H, hH;

		// Weights
		Array<double> wx(ndivx);  wx = 1.0;
		Array<double> wy(ndivy);  wy = 1.0;

		// Blocks
		Array<Mesh::Block> blocks;  blocks.Resize(1);
		blocks[0].Set (c, wx, wy);

		// Generate
		std::clock_t start = std::clock(); // Initial time
		Mesh::Structured m;
		size_t ne = m.Generate (blocks);
		std::clock_t total = std::clock() - start; // Time elapsed
		std::cout << "2D("<<ne<<" elements): Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

		// Output
		m.WriteVTU ("tmesh01_2D.vtu");
		cout << "[1;34mFile <tmesh01_2D.vtu> created[0m" << endl << endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////// 3D /////
	
	{
		// Constants
		double L     = 1.0;   // length
		double H     = 1.0;   // height
		double D     = 1.0;   // depth
		double hL    = L/2.0; // half-length
		double hH    = H/2.0; // half-height
		double hD    = D/2.0; // half-depth
		size_t ndivx = 18;    // divisions along x
		size_t ndivy = 18;    // divisions along y
		size_t ndivz = 18;    // divisions along z

		// Coordinates
		Matrix<double> c(3,20);
		c = 0.,  L,  L, 0.,    0.,  L, L, 0.,    hL,  L, hL, 0.,    hL,  L, hL, 0.,    0.,  L,  L, 0.,
			0., 0.,  H,  H,    0., 0., H,  H,    0., hH,  H, hH,    0., hH,  H, hH,    0., 0.,  H,  H,
		    0., 0., 0., 0.,     D,  D, D,  D,    0., 0., 0., 0.,     D,  D,  D,  D,    hD, hD, hD, hD;

		// Weights
		Array<double> wx(ndivx);  wx = 1.0;
		Array<double> wy(ndivy);  wy = 1.0;
		Array<double> wz(ndivz);  wz = 1.0;

		// Blocks
		Array<Mesh::Block> blocks;  blocks.Resize(1);
		blocks[0].Set (c, wx, wy, wz);

		// Generate
		std::clock_t start = std::clock(); // Initial time
		Mesh::Structured m(/*Is3D*/true);
		size_t ne = m.Generate (blocks);
		std::clock_t total = std::clock() - start; // Time elapsed
		std::cout << "3D:("<<ne<<" elements) Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

		// Output
		m.WriteVTU ("tmesh01_3D.vtu");
		cout << "[1;34mFile <tmesh01_3D.vtu> created[0m" << endl;
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
